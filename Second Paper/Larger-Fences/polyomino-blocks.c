/*
 * polyomino-blocks.c
 *
 * Same sampler and output as polyomino-grid.c, but the grid is stored as
 * 4x4 blocks: each block is a 16-bit pattern of occupied cells (bit
 * 4*row + col), and the blocks of the surrounding box form a 2-dimensional
 * array of uint16_t. The perimeter cells get a second block array.
 *
 * For each of the 2^16 patterns a table (BlockInfo, 32 bytes) is built once
 * at start-up:
 *   - the 4-connected components inside the block, numbered so that the
 *     components touching the block boundary come first (at most 6 of
 *     them, since there are 12 boundary cells),
 *   - a 4-bit label per cell (0xF for empty),
 *   - for each component, which cells of each of the 4 sides it owns
 *     (4 sides x 4 bits),
 *   - which cells of each side are occupied at all,
 *   - whether the pattern is feasible: no component cut off from the block
 *     boundary. For n >= 5 such a component is a whole piece of at most 4
 *     cells, so the polyomino would be disconnected.
 *
 * Connectivity check after a move: if a changed block is infeasible, reject
 * at once. Otherwise run the simultaneous BFS on the graph whose nodes are
 * (block, component) pairs; two nodes in neighbouring blocks are joined when
 * the facing side masks overlap. The BFS never looks inside a block. If all
 * tiles next to x lie in the same component of their block, no search is
 * needed at all.
 *
 * Build:  cc -O2 -o polyomino-blocks polyomino-blocks.c -lm
 * Usage:  ./polyomino-blocks n id [seed [samples]]
 *         Runs until Ctrl+C unless a sample count is given. With the same
 *         seed it writes exactly the same rows as polyomino-grid.
 *         Build with -DCHECK to verify every step against brute force.
 */

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define BURN_IN_FACTOR 4
#define MARGIN 5   /* every tile stays at least MARGIN cells from the edge,
                      so tiles never sit in an edge block */
#ifndef EXTRA_WIDTH /* grid width is about n + 2*MARGIN + EXTRA_WIDTH */
#define EXTRA_WIDTH (n / 2 + 8)
#endif

/* Neighbour order as in the Julia code: +x, -y, -x, +y */
static const int DX[4] = {1, 0, -1, 0};
static const int DY[4] = {0, -1, 0, 1};

static volatile sig_atomic_t stop_requested = 0;
static void on_sigint(int sig) { (void)sig; stop_requested = 1; }

static void *xmalloc(size_t size)
{
    void *p = calloc(1, size);
    if (!p) { fprintf(stderr, "Error: out of memory\n"); exit(1); }
    return p;
}

/* ------------- random numbers: xoshiro256** seeded by splitmix64 --------- */

static uint64_t rs[4];

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t rng_next(void)
{
    uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17;
    rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3];
    rs[2] ^= t; rs[3] = rotl(rs[3], 45);
    return r;
}

static void rng_seed(uint64_t s)
{
    for (int i = 0; i < 4; i++) {
        uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rs[i] = z ^ (z >> 31);
    }
}

/* Uniform integer in [0, k) for k > 0 (Lemire's method, unbiased). */
static inline uint32_t rand_below(uint32_t k)
{
    uint64_t m = (rng_next() >> 32) * (uint64_t)k;
    if ((uint32_t)m < k) {
        uint32_t t = (uint32_t)(-k) % k;
        while ((uint32_t)m < t) m = (rng_next() >> 32) * (uint64_t)k;
    }
    return (uint32_t)(m >> 32);
}

/* ---------------------- per-pattern block tables -------------------------- */

typedef struct {
    uint64_t label;     /* 4 bits per cell: component number, 0xF if empty */
    uint16_t side[8];   /* per component: bit 4*s + i if it owns cell i of side s */
    uint16_t occ_side;  /* bit 4*s + i if cell i of side s is occupied */
    uint8_t nbc;        /* components touching the boundary (numbered first) */
    uint8_t ncomp;      /* all components (at most 8, e.g. a checkerboard) */
    uint8_t feasible;   /* 1 if every component touches the boundary */
    uint8_t pad[3];
} BlockInfo;

_Static_assert(sizeof(BlockInfo) == 32, "BlockInfo should be 32 bytes");

static BlockInfo INFO[65536];

/* Sides in neighbour order (+x, -y, -x, +y) = (E, S, W, N). Cell i of a side
   is counted along the side, so cell i of side s faces cell i of side s^2
   (the opposite side) in the neighbouring block. */
static const int SIDE_CELL[4][4] = {
    {3, 7, 11, 15},    /* E: column 3 */
    {0, 1, 2, 3},      /* S: row 0 */
    {0, 4, 8, 12},     /* W: column 0 */
    {12, 13, 14, 15},  /* N: row 3 */
};

static void build_tables(void)
{
    static const int ring[12] = {0, 1, 2, 3, 7, 11, 15, 14, 13, 12, 8, 4};
    static const int inner[4] = {5, 6, 9, 10};
    for (int P = 0; P < 65536; P++) {
        BlockInfo *bi = &INFO[P];
        int lab[16], nc = 0;
        for (int i = 0; i < 16; i++) lab[i] = -1;
        for (int pass = 0; pass < 2; pass++) {
            const int *cells = pass ? inner : ring;
            for (int s = 0; s < (pass ? 4 : 12); s++) {
                int start = cells[s];
                if (!((P >> start) & 1) || lab[start] >= 0) continue;
                int stack[16], sp = 0;
                stack[sp++] = start; lab[start] = nc;
                while (sp) {
                    int t = stack[--sp], r = t >> 2, col = t & 3;
                    int nbs[4] = {col < 3 ? t + 1 : -1, r > 0 ? t - 4 : -1,
                                  col > 0 ? t - 1 : -1, r < 3 ? t + 4 : -1};
                    for (int d = 0; d < 4; d++)
                        if (nbs[d] >= 0 && ((P >> nbs[d]) & 1) && lab[nbs[d]] < 0) {
                            lab[nbs[d]] = nc; stack[sp++] = nbs[d];
                        }
                }
                nc++;
            }
            if (pass == 0) bi->nbc = (uint8_t)nc;
        }
        bi->ncomp = (uint8_t)nc;
        bi->feasible = nc == bi->nbc;
        bi->label = 0;
        for (int i = 0; i < 16; i++)
            bi->label |= (uint64_t)(lab[i] < 0 ? 15 : lab[i]) << (4 * i);
        for (int s = 0; s < 4; s++)
            for (int i = 0; i < 4; i++) {
                int cell = SIDE_CELL[s][i];
                if (lab[cell] < 0) continue;
                bi->side[lab[cell]] |= (uint16_t)(1u << (4 * s + i));
                bi->occ_side |= (uint16_t)(1u << (4 * s + i));
            }
    }
}

/* ------------------------------- chain state ------------------------------ */

typedef struct { int32_t x, y; } Cell;

typedef struct {
    int n, W, NB;            /* polyomino size, grid W x W cells = NB x NB blocks */
    int boff[4];             /* block index offsets of the 4 neighbouring blocks */
    uint16_t *tblk, *pblk;   /* tile / perimeter pattern of each block */
    int32_t *pos;            /* index of the cell in tiles[] or peri[], by y*W + x */
    Cell *tiles; int ntiles;
    Cell *peri;  int nperi;
    /* BFS: visit[b] = gen << 32 | 3-bit mark per component of block b */
    uint64_t *visit; uint32_t gen;
    int32_t *queue[4];
} Chain;

static inline int bidx(const Chain *c, int x, int y) { return (y >> 2) * c->NB + (x >> 2); }
static inline int bbit(int x, int y) { return ((y & 3) << 2) | (x & 3); }
static inline int is_tile(const Chain *c, int x, int y) { return (c->tblk[bidx(c, x, y)] >> bbit(x, y)) & 1; }
static inline int is_peri(const Chain *c, int x, int y) { return (c->pblk[bidx(c, x, y)] >> bbit(x, y)) & 1; }
static inline int same_cell(Cell a, Cell b) { return a.x == b.x && a.y == b.y; }

static inline void tile_add(Chain *c, Cell p)
{
    c->tblk[bidx(c, p.x, p.y)] |= (uint16_t)(1u << bbit(p.x, p.y));
    c->pos[p.y * c->W + p.x] = c->ntiles; c->tiles[c->ntiles++] = p;
}

static inline void tile_remove(Chain *c, Cell p)
{
    int k = c->pos[p.y * c->W + p.x];
    Cell last = c->tiles[--c->ntiles];
    c->tiles[k] = last; c->pos[last.y * c->W + last.x] = k;
    c->tblk[bidx(c, p.x, p.y)] &= (uint16_t)~(1u << bbit(p.x, p.y));
}

static inline void peri_add(Chain *c, Cell p)
{
    c->pblk[bidx(c, p.x, p.y)] |= (uint16_t)(1u << bbit(p.x, p.y));
    c->pos[p.y * c->W + p.x] = c->nperi; c->peri[c->nperi++] = p;
}

static inline void peri_remove(Chain *c, Cell p)
{
    int k = c->pos[p.y * c->W + p.x];
    Cell last = c->peri[--c->nperi];
    c->peri[k] = last; c->pos[last.y * c->W + last.x] = k;
    c->pblk[bidx(c, p.x, p.y)] &= (uint16_t)~(1u << bbit(p.x, p.y));
}

static Chain *chain_new(int n)
{
    Chain *c = xmalloc(sizeof *c);
    c->n = n;
    c->NB = (n + 2 * MARGIN + EXTRA_WIDTH + 3) / 4;
    c->W = 4 * c->NB;
    int boff[4] = {1, -c->NB, -1, c->NB};
    memcpy(c->boff, boff, sizeof boff);
    size_t blocks = (size_t)c->NB * c->NB;
    c->tblk  = xmalloc(blocks * sizeof(uint16_t));
    c->pblk  = xmalloc(blocks * sizeof(uint16_t));
    c->visit = xmalloc(blocks * sizeof(uint64_t));
    c->pos   = xmalloc((size_t)c->W * c->W * sizeof(int32_t));
    c->tiles = xmalloc((size_t)(n + 4) * sizeof(Cell));
    c->peri  = xmalloc((size_t)(4 * n + 8) * sizeof(Cell));
    for (int l = 0; l < 4; l++) c->queue[l] = xmalloc((size_t)(n + 8) * sizeof(int32_t));
    return c;
}

static void clear_blocks(Chain *c)
{
    for (int k = 0; k < c->ntiles; k++) c->tblk[bidx(c, c->tiles[k].x, c->tiles[k].y)] = 0;
    for (int k = 0; k < c->nperi; k++)  c->pblk[bidx(c, c->peri[k].x, c->peri[k].y)] = 0;
}

/* Straight horizontal line of n tiles in the middle of the grid. */
static void chain_init(Chain *c)
{
    clear_blocks(c);
    c->ntiles = c->nperi = 0;
    int x0 = (c->W - c->n) / 2, y0 = c->W / 2;
    for (int i = 0; i < c->n; i++) tile_add(c, (Cell){x0 + i, y0});
    for (int k = 0; k < c->ntiles; k++)
        for (int d = 0; d < 4; d++) {
            Cell t = {c->tiles[k].x + DX[d], c->tiles[k].y + DY[d]};
            if (!is_tile(c, t.x, t.y) && !is_peri(c, t.x, t.y)) peri_add(c, t);
        }
}

static long long recentre_count = 0;

/* Move the polyomino (and its perimeter) to the middle of the grid. */
static void recentre(Chain *c)
{
    recentre_count++;
    int W = c->W, minx = W, maxx = -1, miny = W, maxy = -1;
    for (int k = 0; k < c->ntiles; k++) {
        Cell p = c->tiles[k];
        if (p.x < minx) minx = p.x;
        if (p.x > maxx) maxx = p.x;
        if (p.y < miny) miny = p.y;
        if (p.y > maxy) maxy = p.y;
    }
    int sx = (W - (maxx - minx + 1)) / 2 - minx;
    int sy = (W - (maxy - miny + 1)) / 2 - miny;
    clear_blocks(c);
    for (int k = 0; k < c->ntiles; k++) {
        Cell *p = &c->tiles[k];
        p->x += sx; p->y += sy;
        c->tblk[bidx(c, p->x, p->y)] |= (uint16_t)(1u << bbit(p->x, p->y));
        c->pos[p->y * W + p->x] = k;
    }
    for (int k = 0; k < c->nperi; k++) {
        Cell *p = &c->peri[k];
        p->x += sx; p->y += sy;
        c->pblk[bidx(c, p->x, p->y)] |= (uint16_t)(1u << bbit(p->x, p->y));
        c->pos[p->y * W + p->x] = k;
    }
}

static inline int uf_find(const int *parent, int a)
{
    while (parent[a] != a) a = parent[a];
    return a;
}

/*
 * Simultaneous BFS on (block, component) nodes, starting from the nodes of
 * the k (2..4) tiles in start[]. Nodes are encoded as block * 8 + component.
 * Every reached node is marked with the label of the start it was reached
 * from; meeting a node with a mark from another group merges the groups.
 * One node is taken from each label's queue per round. Returns 1 as soon as
 * all groups have merged, 0 as soon as some group has run out of nodes while
 * others remain.
 */
static int bfs_blocks(Chain *c, const Cell *start, int k)
{
    if (++c->gen == 0) {
        memset(c->visit, 0, (size_t)c->NB * c->NB * sizeof(uint64_t));
        c->gen = 1;
    }
    const uint64_t tag = (uint64_t)c->gen << 32, TAGMASK = 0xFFFFFFFF00000000ULL;
    uint64_t *visit = c->visit;
    const uint16_t *tblk = c->tblk;

    int node[4], ns = 0;
    for (int i = 0; i < k; i++) {
        int b = bidx(c, start[i].x, start[i].y);
        int comp = (int)(INFO[tblk[b]].label >> (4 * bbit(start[i].x, start[i].y))) & 15;
        int nd = b * 8 + comp, dup = 0;
        for (int j = 0; j < ns; j++) dup |= node[j] == nd;
        if (!dup) node[ns++] = nd;
    }
    if (ns <= 1) return 1;   /* all joined inside one block */

    int head[4], tail[4], parent[4], groups = ns;
    for (int l = 0; l < ns; l++) {
        int b = node[l] >> 3, comp = node[l] & 7;
        c->queue[l][0] = node[l]; head[l] = 0; tail[l] = 1; parent[l] = l;
        uint64_t v = visit[b];
        if ((v & TAGMASK) != tag) v = tag;
        visit[b] = v | ((uint64_t)(l + 1) << (3 * comp));
    }
    for (;;) {
        for (int l = 0; l < ns; l++) {
            if (head[l] == tail[l]) continue;
            int nd = c->queue[l][head[l]++], b = nd >> 3, comp = nd & 7;
            unsigned sm = INFO[tblk[b]].side[comp];
            for (int s = 0; s < 4; s++) {
                unsigned m = (sm >> (4 * s)) & 15;
                if (!m) continue;
                int b2 = b + c->boff[s], so = s ^ 2;
                const BlockInfo *bj = &INFO[tblk[b2]];
                unsigned f = m & (bj->occ_side >> (4 * so));
                while (f) {
                    int i = __builtin_ctz(f);
                    f &= f - 1;
                    int c2 = (int)(bj->label >> (4 * SIDE_CELL[so][i])) & 15;
                    uint64_t v = visit[b2];
                    if ((v & TAGMASK) != tag) v = tag;
                    int mk = (int)(v >> (3 * c2)) & 7;
                    if (!mk) {
                        visit[b2] = v | ((uint64_t)(l + 1) << (3 * c2));
                        c->queue[l][tail[l]++] = b2 * 8 + c2;
                    } else {
                        int ra = uf_find(parent, l), rb = uf_find(parent, mk - 1);
                        if (ra != rb) {
                            parent[rb] = ra;
                            if (--groups == 1) return 1;
                        }
                    }
                }
            }
        }
        /* a group whose queues are all empty has explored its whole
           component without meeting the others: disconnected */
        for (int r = 0; r < ns; r++) {
            if (parent[r] != r) continue;
            int alive = 0;
            for (int l = 0; l < ns && !alive; l++)
                alive = head[l] < tail[l] && uf_find(parent, l) == r;
            if (!alive) return 0;
        }
    }
}

#ifdef CHECK
/* Reference: plain cell flood fill from start[0], all others must be reached. */
static int check_connected(Chain *c, const Cell *start, int k)
{
    uint8_t *seen = xmalloc((size_t)c->W * c->W);
    Cell *stack = xmalloc((size_t)(c->n + 4) * sizeof(Cell));
    int sp = 0, ok = 1, W = c->W;
    stack[sp++] = start[0]; seen[start[0].y * W + start[0].x] = 1;
    while (sp) {
        Cell t = stack[--sp];
        for (int d = 0; d < 4; d++) {
            Cell u = {t.x + DX[d], t.y + DY[d]};
            if (is_tile(c, u.x, u.y) && !seen[u.y * W + u.x]) { seen[u.y * W + u.x] = 1; stack[sp++] = u; }
        }
    }
    for (int i = 1; i < k; i++) ok &= seen[start[i].y * W + start[i].x];
    free(seen); free(stack);
    return ok;
}

/* Arrays, pos[] and block bits agree, perimeter is exact. */
static void check_state(Chain *c)
{
    int W = c->W, nt = 0, np = 0;
    for (int k = 0; k < c->ntiles; k++) {
        Cell p = c->tiles[k];
        if (!is_tile(c, p.x, p.y) || is_peri(c, p.x, p.y) || c->pos[p.y * W + p.x] != k) goto bad;
        if (p.x < MARGIN - 1 || p.x > W - MARGIN || p.y < MARGIN - 1 || p.y > W - MARGIN) goto bad;
    }
    for (int k = 0; k < c->nperi; k++) {
        Cell p = c->peri[k];
        if (!is_peri(c, p.x, p.y) || c->pos[p.y * W + p.x] != k) goto bad;
    }
    for (int y = 1; y < W - 1; y++)
        for (int x = 1; x < W - 1; x++) {
            int adj = 0;
            for (int d = 0; d < 4; d++) adj |= is_tile(c, x + DX[d], y + DY[d]);
            if (is_tile(c, x, y)) nt++;
            else if (adj != is_peri(c, x, y)) goto bad;
            if (is_peri(c, x, y)) np++;
        }
    if (nt != c->ntiles || np != c->nperi || nt != c->n) goto bad;
    return;
bad:
    fprintf(stderr, "CHECK: inconsistent chain state\n");
    exit(2);
}

/* Block tables against a direct recount for a pattern. */
static void check_tables(void)
{
    for (int P = 0; P < 65536; P++) {
        const BlockInfo *bi = &INFO[P];
        unsigned all = 0, occ = 0;
        for (int c = 0; c < bi->ncomp; c++) {
            if ((bi->side[c] & all) || (c < bi->nbc) != (bi->side[c] != 0)) goto bad;
            all |= bi->side[c];
        }
        for (int s = 0; s < 4; s++)
            for (int i = 0; i < 4; i++)
                if ((P >> SIDE_CELL[s][i]) & 1) occ |= 1u << (4 * s + i);
        if (all != occ || occ != bi->occ_side || bi->nbc > 6 || bi->ncomp > 8) goto bad;
        for (int i = 0; i < 16; i++) {
            int lab = (int)(bi->label >> (4 * i)) & 15;
            if (((P >> i) & 1) != (lab != 15)) goto bad;
            /* neighbours inside the block share the label */
            if (lab != 15 && (i & 3) < 3 && ((P >> (i + 1)) & 1) && lab != (int)((bi->label >> (4 * i + 4)) & 15)) goto bad;
            if (lab != 15 && i < 12 && ((P >> (i + 4)) & 1) && lab != (int)((bi->label >> (4 * i + 16)) & 15)) goto bad;
        }
    }
    return;
bad:
    fprintf(stderr, "CHECK: block tables inconsistent\n");
    exit(2);
}
#endif

/*
 * One Metropolis step (same moves as the Julia code): remove a random tile
 * x, add a random perimeter cell y, keep the move if the result is
 * connected, otherwise roll it back. Returns 1 if accepted.
 */
static int metropolis_step(Chain *c)
{
    /* Step 1: remove a random tile x */
    Cell x = c->tiles[rand_below((uint32_t)c->ntiles)];
    tile_remove(c, x);
    peri_add(c, x);

    Cell del[4]; int ndel = 0;
    for (int d = 0; d < 4; d++) {
        Cell t = {x.x + DX[d], x.y + DY[d]};
        if (!is_peri(c, t.x, t.y)) continue;
        int adj = 0;
        for (int e = 0; e < 4; e++) adj |= is_tile(c, t.x + DX[e], t.y + DY[e]);
        if (!adj) del[ndel++] = t;
    }
    for (int i = 0; i < ndel; i++) peri_remove(c, del[i]);

    /* Step 2: add a random perimeter cell y */
    Cell y = c->peri[rand_below((uint32_t)c->nperi)];
    peri_remove(c, y);
    tile_add(c, y);

    Cell add[4]; int nadd = 0;
    for (int d = 0; d < 4; d++) {
        Cell t = {y.x + DX[d], y.y + DY[d]};
        if (!is_tile(c, t.x, t.y) && !is_peri(c, t.x, t.y)) add[nadd++] = t;
    }
    for (int i = 0; i < nadd; i++) peri_add(c, add[i]);

    /* Step 3: connectivity check. y == x puts the old polyomino back. */
    int connected = 1;
    if (!same_cell(x, y)) {
        Cell nb[4]; int k = 0;
        for (int d = 0; d < 4; d++)
            if (is_tile(c, x.x + DX[d], x.y + DY[d])) nb[k++] = (Cell){x.x + DX[d], x.y + DY[d]};
        if (c->n >= 5 && (!INFO[c->tblk[bidx(c, x.x, x.y)]].feasible ||
                          !INFO[c->tblk[bidx(c, y.x, y.y)]].feasible))
            connected = 0;   /* a piece of <= 4 cells is cut off inside a block */
        else if (k > 1)
            connected = bfs_blocks(c, nb, k);
#ifdef CHECK
        if (connected != (k > 1 ? check_connected(c, nb, k) : 1)) {
            fprintf(stderr, "CHECK: connectivity mismatch\n");
            exit(2);
        }
#endif
    }

    if (connected) {
        if (y.x < MARGIN || y.x >= c->W - MARGIN || y.y < MARGIN || y.y >= c->W - MARGIN)
            recentre(c);
        return 1;
    }

    /* Rollback (same order of array operations as the Julia code) */
    tile_remove(c, y);
    for (int i = 0; i < nadd; i++) peri_remove(c, add[i]);
    peri_add(c, y);
    for (int i = 0; i < ndel; i++) peri_add(c, del[i]);
    peri_remove(c, x);
    tile_add(c, x);
    return 0;
}

/* Advance by `steps` accepted or rejected steps; returns 0 if interrupted. */
static int advance_chain(Chain *c, long long steps)
{
    for (long long done = 0; done < steps; ) {
        if (stop_requested) return 0;
	metropolis_step(c);
        done += 1;
#ifdef CHECK
        check_state(c);
#endif
    }
    return 1;
}

/* --------------- diameter: convex hull + rotating calipers --------------- */

typedef struct { int x, y, i; } Pt;   /* doubled coordinates, original index */

static int pt_cmp(const void *a, const void *b)
{
    const Pt *p = a, *q = b;
    if (p->x != q->x) return p->x < q->x ? -1 : 1;
    if (p->y != q->y) return p->y < q->y ? -1 : 1;
    return 0;
}

static inline long long cross(Pt o, Pt a, Pt b)
{
    return (long long)(a.x - o.x) * (b.y - o.y) - (long long)(a.y - o.y) * (b.x - o.x);
}

static inline long long dist2(Pt a, Pt b)
{
    long long dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

/* Andrew's monotone chain. Sorts pts in place, writes the hull in
   counter-clockwise order to hull[] (room for m+1 points), drops collinear
   points (they can never be an endpoint of a farthest pair). */
static int convex_hull(Pt *pts, int m, Pt *hull)
{
    qsort(pts, (size_t)m, sizeof(Pt), pt_cmp);
    if (m <= 2) { memcpy(hull, pts, (size_t)m * sizeof(Pt)); return m; }
    int h = 0;
    for (int i = 0; i < m; i++) {
        while (h >= 2 && cross(hull[h - 2], hull[h - 1], pts[i]) <= 0) h--;
        hull[h++] = pts[i];
    }
    for (int i = m - 2, lo = h + 1; i >= 0; i--) {
        while (h >= lo && cross(hull[h - 2], hull[h - 1], pts[i]) <= 0) h--;
        hull[h++] = pts[i];
    }
    return h - 1;
}

typedef struct { int a, b; } Pair;

static int pair_cmp(const void *p, const void *q)
{
    const Pair *u = p, *v = q;
    if (u->a != v->a) return u->a < v->a ? -1 : 1;
    if (u->b != v->b) return u->b < v->b ? -1 : 1;
    return 0;
}

static long long best_d;
static int npairs;

static inline void consider(Pt p, Pt q, Pair *pairs)
{
    long long d = dist2(p, q);
    if (d < best_d) return;
    if (d > best_d) { best_d = d; npairs = 0; }
    pairs[npairs].a = p.i < q.i ? p.i : q.i;
    pairs[npairs].b = p.i < q.i ? q.i : p.i;
    npairs++;
}

/* All pairs (original indices, a < b) of hull points at maximal distance,
   sorted and without duplicates. Rotating calipers visits every antipodal
   pair, including both ends of an edge parallel to the current one, so no
   tie is lost. pairs[] needs room for 4h + 4 entries. */
static int diametral_pairs(const Pt *H, int h, Pair *pairs)
{
    best_d = -1; npairs = 0;
    if (h < 2) return 0;
    if (h == 2) { consider(H[0], H[1], pairs); return npairs; }
    int j = 1;
    for (int i = 0; i < h; i++) {
        int i2 = (i + 1) % h;
        while (cross(H[i], H[i2], H[(j + 1) % h]) > cross(H[i], H[i2], H[j])) j = (j + 1) % h;
        consider(H[i], H[j], pairs); consider(H[i2], H[j], pairs);
        int jn = (j + 1) % h;
        if (cross(H[i], H[i2], H[jn]) == cross(H[i], H[i2], H[j])) {
            consider(H[i], H[jn], pairs); consider(H[i2], H[jn], pairs);
        }
    }
    qsort(pairs, (size_t)npairs, sizeof(Pair), pair_cmp);
    int u = 0;
    for (int k = 0; k < npairs; k++)
        if (u == 0 || pair_cmp(&pairs[k], &pairs[u - 1]) != 0) pairs[u++] = pairs[k];
    return npairs = u;
}

typedef struct { Pt *pts, *orig, *hull; Pair *pairs; } Scratch;

static Scratch scratch_new(int n)
{
    Scratch s;
    int m = 4 * n + 4;
    s.pts = xmalloc((size_t)m * sizeof(Pt));
    s.orig = xmalloc((size_t)m * sizeof(Pt));
    s.hull = xmalloc((size_t)(m + 1) * sizeof(Pt));
    s.pairs = xmalloc((size_t)(4 * m + 8) * sizeof(Pair));
    return s;
}

/* Write one row per diametral pair: n, folded x, folded y, draw_size. */
static void write_sample(FILE *out, Chain *c, Scratch *s)
{
    int m = 0;
    for (int k = 0; k < c->ntiles; k++) {
        Cell p = c->tiles[k];
        for (int d = 0; d < 4; d++)
            if (!is_tile(c, p.x + DX[d], p.y + DY[d])) {
                Pt q = {2 * p.x + DX[d], 2 * p.y + DY[d], m};
                s->orig[m] = q; s->pts[m++] = q;
            }
    }
    int h = convex_hull(s->pts, m, s->hull);
    int np = diametral_pairs(s->hull, h, s->pairs);
#ifdef CHECK
    long long best = -1; int cnt = 0;
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++) {
            long long d = dist2(s->orig[i], s->orig[j]);
            if (d > best) { best = d; cnt = 0; }
            if (d == best) cnt++;
        }
    if (cnt != np || best != best_d) { fprintf(stderr, "CHECK: diameter mismatch\n"); exit(2); }
#endif
    for (int k = 0; k < np; k++) {
        int ax = abs(s->orig[s->pairs[k].b].x - s->orig[s->pairs[k].a].x);
        int ay = abs(s->orig[s->pairs[k].b].y - s->orig[s->pairs[k].a].y);
        if (ay > ax) { int t = ax; ax = ay; ay = t; }
        /* coordinates are doubled; a diametral vector always joins two
           midpoints of the same kind, so both components are even */
        if ((ax | ay) & 1) {
            fprintf(stderr, "Error: diametral vector with a half-integer component\n");
            exit(1);
        }
        fprintf(out, "%d,%d,%d,%d\n", c->n, ax / 2, ay / 2, np);
    }
}

/* ----------------------------------- main ---------------------------------- */

int main(int argc, char **argv)
{
    if (argc < 3 || argc > 5) {
        printf("Usage: %s n id [seed [samples]]\n", argv[0]);
        return 1;
    }
    char *end;
    long n = strtol(argv[1], &end, 10);
    if (*end || end == argv[1]) { printf("Error: n must be an integer.\n"); return 1; }
    if (n < 1 || n > 100000) { printf("Error: n must be between 1 and 100000 (got %ld).\n", n); return 1; }
    const char *id = argv[2];
    if (!*id || strchr(id, '/') || strchr(id, '\\')) {
        printf("Error: id must be non-empty and must not contain '/' or '\\'.\n");
        return 1;
    }
    uint64_t seed = argc > 3 ? strtoull(argv[3], NULL, 10)
                             : (uint64_t)time(NULL) ^ ((uint64_t)getpid() << 32);
    long long max_samples = argc > 4 ? strtoll(argv[4], NULL, 10) : 0;
    rng_seed(seed);

    if (mkdir("data", 0777) != 0 && errno != EEXIST) {
        fprintf(stderr, "Error: could not create data directory: %s\n", strerror(errno));
        return 1;
    }
    char path[4096];
    snprintf(path, sizeof path, "data/vector-data-%s.csv", id);
    int write_header = access(path, F_OK) != 0;
    FILE *out = fopen(path, "a");
    if (!out) { fprintf(stderr, "Error: cannot open %s: %s\n", path, strerror(errno)); return 1; }
    if (write_header) fprintf(out, "n,x,y,draw_size\n");

    signal(SIGINT, on_sigint);
    printf("n=%ld  id=%s  seed=%llu  burn-in=%lld per sample  writing to %s (Ctrl+C to stop)\n",
           n, id, (unsigned long long)seed, (long long)BURN_IN_FACTOR * n * n, path);
    fflush(stdout);

    build_tables();
#ifdef CHECK
    check_tables();
#endif
    Chain *c = chain_new((int)n);
    Scratch s = scratch_new((int)n);
    long long sampled = 0;
    chain_init(c);
    while (!max_samples || sampled < max_samples) {
        if (!advance_chain(c, (long long)BURN_IN_FACTOR * n * n)) break;
        write_sample(out, c, &s);
        fflush(out);
        if (++sampled % 100 == 0) { printf("n=%ld  id=%s  samples=%lld\n", n, id, sampled); fflush(stdout); }
    }
    fclose(out);
    if (stop_requested) printf("\nStopped by user (n=%ld, id=%s, samples=%lld).\n", n, id, sampled);
#ifdef CHECK
    fprintf(stderr, "CHECK: all checks passed (%lld samples, %lld recentrings)\n", sampled, recentre_count);
#endif
    return 0;
}
