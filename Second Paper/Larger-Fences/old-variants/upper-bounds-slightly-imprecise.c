/* Fmax + Dmax bound U(n) and the linearised bound S3 for n-omino fences.
   Re-implementation of upper-bounds.py (compute_ref_vectors) in C,
   provided by Claude from Anthropic.
   Usage: upper-bounds n target
   Pass 1: u_max, K_max (sum of max corrections) -> U(n), A0, B0.
   Pass 2: S3 = sum_i max_j (A0 |v_ij| + Delta_ij) + B0.
   All coordinates doubled; corrections handled as integers in units of 1/8. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int n;
static int GW, GH;                 /* Redelmeier grid */
static unsigned char *seen;
static int poly[64];
static long long nfixed = 0, nfree = 0;
static int pass;
static long double A0, B0;

/* pass 1 accumulators */
static long long *len2count;       /* histogram of longest len2 (doubled coords) */
static int maxlen2_global;
static long long sumK8 = 0;        /* sum of max corrections *8 */
static long long sumKtop8 = 0;     /* sum of max correction among longest vectors *8 */
/* pass 2 accumulators */
static long double s3sum = 0, s3c = 0;   /* Kahan */
static long long nonlongest = 0;
static long double maxgainshort = -1e30L; /* max over pieces of best - top */

static void kahan(long double v) {
    long double y = v - s3c, t = s3sum + y;
    s3c = (t - s3sum) - y; s3sum = t;
}

/* ---- boundary and reference vectors ---- */
static int cx[64], cy[64];

typedef struct { int x, y; } V;

static void process_path(V *p, int L, int *maxlen2, int *maxD8, int *topD8, long double *best)
{
    static long long au[256], pre[512];
    for (int k = 0; k < L; k++) {
        V a = p[k], b = p[(k + 1) % L];
        au[k] = (long long)(b.x - a.x) * (a.y + b.y);
    }
    pre[0] = 0;
    for (int k = 0; k < 2 * L; k++) pre[k + 1] = pre[k] + au[k % L];
    for (int i = 0; i < L; i++) {
        V a = p[i], b = p[(i + 1) % L];
        int xS2 = a.x + b.x, yS2 = a.y + b.y;   /* midpoint *2 in doubled coords */
        for (int j = i + 1; j < i + L; j++) {
            V c = p[j % L], d = p[(j + 1) % L];
            int xE2 = c.x + d.x, yE2 = c.y + d.y;
            /* along (doubled coords) = (au_i + au_j)/2 + pre[j] - pre[i+1] */
            long long along2 = au[i] + au[j % L] + 2 * (pre[j] - pre[i + 1]);  /* *2 */
            long long dx2 = xE2 - xS2, dy2 = yE2 - yS2;           /* *2 */
            long long aref4 = dx2 * (yE2 + yS2);                    /* *4 of doubled */
            /* P8 = along_d - aref_d ; along_d = along2/2, aref_d = aref4/4 */
            long long num = 2 * along2 - aref4;                     /* = 4*P8 */
            if (num % 4) { fprintf(stderr, "non-integer P8\n"); exit(1); }
            long long P8 = num / 4;
            int D8 = (int)(-P8);
            long long l2 = dx2 * dx2 + dy2 * dy2;                   /* = 4*len2_d */
            if (l2 % 4) { fprintf(stderr, "odd l2\n"); exit(1); }
            int len2 = (int)(l2 / 4);                               /* doubled-coord len^2 = 4|v|^2 */
            if (pass == 1) {
                if (len2 > *maxlen2) { *maxlen2 = len2; *topD8 = D8; }
                else if (len2 == *maxlen2 && D8 > *topD8) *topD8 = D8;
                if (D8 > *maxD8) *maxD8 = D8;
            } else {
                long double val = A0 * sqrtl((long double)len2) / 2 + D8 / 8.0L;
                if (val > *best) *best = val;
                if (len2 > *maxlen2) { *maxlen2 = len2; *topD8 = D8; }
                else if (len2 == *maxlen2 && D8 > *topD8) *topD8 = D8;
            }
        }
    }
}

static void process(void)
{
    /* coordinates */
    int W = GW;
    int minx = 1 << 30, miny = 1 << 30, maxx = -(1 << 30), maxy = -(1 << 30);
    for (int k = 0; k < n; k++) {
        cx[k] = poly[k] % W; cy[k] = poly[k] / W;
        if (cx[k] < minx) minx = cx[k];
        if (cy[k] < miny) miny = cy[k];
        if (cx[k] > maxx) maxx = cx[k];
        if (cy[k] > maxy) maxy = cy[k];
    }
    nfixed++;
    /* canonical check: identity code vs 8 transforms */
    int code0[64], code[64];
    for (int k = 0; k < n; k++) code0[k] = (cx[k] - minx) * 64 + (cy[k] - miny);
    for (int a = 1; a < n; a++) { int t = code0[a], b = a - 1; while (b >= 0 && code0[b] > t) { code0[b + 1] = code0[b]; b--; } code0[b + 1] = t; }
    for (int s = 1; s < 8; s++) {
        int tx[64], ty[64], mx = 1 << 30, my = 1 << 30;
        for (int k = 0; k < n; k++) {
            int x = cx[k], y = cy[k], X, Y;
            switch (s) {
            case 1: X = -y; Y = x; break;
            case 2: X = -x; Y = -y; break;
            case 3: X = y; Y = -x; break;
            case 4: X = y; Y = x; break;
            case 5: X = -x; Y = y; break;
            case 6: X = -y; Y = -x; break;
            default: X = x; Y = -y; break;
            }
            tx[k] = X; ty[k] = Y;
            if (X < mx) mx = X;
            if (Y < my) my = Y;
        }
        for (int k = 0; k < n; k++) code[k] = (tx[k] - mx) * 64 + (ty[k] - my);
        for (int a = 1; a < n; a++) { int t = code[a], b = a - 1; while (b >= 0 && code[b] > t) { code[b + 1] = code[b]; b--; } code[b + 1] = t; }
        int cmp = 0;
        for (int k = 0; k < n && !cmp; k++) cmp = (code[k] < code0[k]) ? -1 : (code[k] > code0[k]) ? 1 : 0;
        if (cmp < 0) return;
    }
    nfree++;

    /* occupancy bitmap with margin 1 */
    int w = maxx - minx + 3, h = maxy - miny + 3;
    static unsigned char occ[64 * 64];
    memset(occ, 0, (size_t)w * h);
    for (int k = 0; k < n; k++) occ[(cy[k] - miny + 1) * w + (cx[k] - minx + 1)] = 1;
    /* vertices in doubled coords: cell (X,Y) local (1..), corners at 2X+-1.
       vertex index: ((vx+1)/2, (vy+1)/2) with vx in [-1, 2w-1] -> (0..w) */
    int vw = w + 1, vh = h + 1;
    static int succ[65 * 65][2];
    static unsigned char nsucc[65 * 65];
    memset(nsucc, 0, (size_t)vw * vh);
    static const int DX[4] = {0, 1, 0, -1}, DY[4] = {-1, 0, 1, 0};
    int startv = -1; int sx = 0, sy = 0;
    for (int k = 0; k < n; k++) {
        int X = cx[k] - minx + 1, Y = cy[k] - miny + 1;
        for (int d = 0; d < 4; d++) {
            int XX = X + DX[d], YY = Y + DY[d];
            if (occ[YY * w + XX]) continue;
            int mx2 = 2 * X + DX[d], my2 = 2 * Y + DY[d];
            int px = -DY[d], py = DX[d];
            int s_x = mx2 + px, s_y = my2 + py, e_x = mx2 - px, e_y = my2 - py;
            int si = ((s_y + 1) / 2) * vw + (s_x + 1) / 2;
            int ei = ((e_y + 1) / 2) * vw + (e_x + 1) / 2;
            succ[si][nsucc[si]++] = ei;
            if (startv < 0 || s_x < sx || (s_x == sx && s_y < sy)) { startv = si; sx = s_x; sy = s_y; }
        }
    }
    if (nsucc[startv] != 1) { fprintf(stderr, "start has %d successors\n", nsucc[startv]); exit(1); }
    V path[256]; int L = 0;
    int prev = startv, cur = succ[startv][0];
    path[L++] = (V){sx, sy};
    while (cur != startv) {
        int cxv = (cur % vw) * 2 - 1, cyv = (cur / vw) * 2 - 1;
        path[L++] = (V){cxv, cyv};
        int nxt;
        if (nsucc[cur] == 1) nxt = succ[cur][0];
        else {
            int pxv = (prev % vw) * 2 - 1, pyv = (prev / vw) * 2 - 1;
            int ddx = cxv - pxv, ddy = cyv - pyv;
            int rdx = ddy, rdy = -ddx;
            int nx = cxv + rdx, ny = cyv + rdy;
            nxt = ((ny + 1) / 2) * vw + (nx + 1) / 2;
            if (succ[cur][0] != nxt && succ[cur][1] != nxt) { fprintf(stderr, "right-hand rule failed\n"); exit(1); }
        }
        prev = cur; cur = nxt;
        if (L > 250) { fprintf(stderr, "path too long\n"); exit(1); }
    }
    /* reflected reversed path */
    V rpath[256];
    for (int k = 0; k < L; k++) { V q = path[L - 1 - k]; rpath[k] = (V){q.y, q.x}; }

    int maxlen2 = -1, maxD8 = -(1 << 30), topD8 = -(1 << 30);
    long double best = -1e30L;
    process_path(path, L, &maxlen2, &maxD8, &topD8, &best);
    process_path(rpath, L, &maxlen2, &maxD8, &topD8, &best);
    if (pass == 1) {
        len2count[maxlen2]++;
        if (maxlen2 > maxlen2_global) maxlen2_global = maxlen2;
        sumK8 += maxD8;
        sumKtop8 += topD8;
    } else {
        kahan(best);
        long double top = A0 * sqrtl((long double)maxlen2) / 2 + topD8 / 8.0L;
        if (best > top + 1e-9L) nonlongest++;
        if (best - top > maxgainshort) maxgainshort = best - top;
    }
}

static void redel(int *untried, int nu, int depth)
{
    int newu[512];
    while (nu > 0) {
        int c = untried[--nu];
        poly[depth] = c;
        if (depth + 1 == n) process();
        else {
            int nu2 = nu;
            memcpy(newu, untried, sizeof(int) * nu);
            int added[4], na = 0;
            int nb[4] = {c + 1, c - 1, c + GW, c - GW};
            for (int k = 0; k < 4; k++)
                if (!seen[nb[k]]) { seen[nb[k]] = 1; newu[nu2++] = nb[k]; added[na++] = nb[k]; }
            redel(newu, nu2, depth + 1);
            for (int k = 0; k < na; k++) seen[added[k]] = 0;
        }
    }
}

static void enumerate(void)
{
    GW = 2 * n + 3; GH = n + 2;
    seen = calloc((size_t)GW * GH, 1);
    /* forbidden: y<0 (row 0 is a guard row), y==0 && x<0; columns 0 and GW-1 guards */
    int oy = 1, ox = n + 1;
    for (int y = 0; y < GH; y++)
        for (int x = 0; x < GW; x++) {
            int X = x - ox, Y = y - oy;
            if (Y < 0 || (Y == 0 && X < 0) || x == 0 || x == GW - 1 || y == GH - 1) seen[y * GW + x] = 1;
        }
    int origin = oy * GW + ox;
    seen[origin] = 1;
    int u[1] = {origin};
    nfixed = nfree = 0;
    redel(u, 1, 0);
    free(seen);
}

int main(int argc, char **argv)
{
    n = atoi(argv[1]);
    long double target = strtold(argv[2], NULL);
    len2count = calloc(100000, sizeof(long long));
    pass = 1;
    enumerate();
    long double k = nfree;
    long double cst = 1.0L / (tanl(acosl(-1.0L) / k) * 4 * k);
    long double u = 0;
    for (int l = 0; l <= maxlen2_global; l++) if (len2count[l]) u += len2count[l] * sqrtl((long double)l) / 2;
    long double F = u * u * cst, K = sumK8 / 8.0L, Ktop = sumKtop8 / 8.0L;
    printf("n=%d fixed=%lld free=%lld\n", n, nfixed, nfree);
    printf("u_max=%.12Lf K_max=%.3Lf F=%.6Lf F+K=%.6Lf U=%.0Lf\n", u, K, F, F + K, floorl(F + K));
    printf("Delta_top=%.3Lf  F+Delta_top=%.6Lf floor=%.0Lf\n", Ktop, F + Ktop, floorl(F + Ktop));
    long double ucrit = sqrtl((target - K) / cst);
    long double Fcrit = ucrit * ucrit * cst;
    A0 = (F - Fcrit) / (u - ucrit);
    B0 = F - u * A0;
    printf("target=%.3Lf u_crit=%.6Lf A0=%.9Lf B0=%.6Lf gap_mid=%.6Lf\n", target, ucrit, A0, B0,
           (Fcrit + F) / 2 - ((ucrit + u) / 2) * ((ucrit + u) / 2) * cst);
    pass = 2;
    enumerate();
    long double S3 = s3sum + B0;
    printf("S3=%.6Lf floor=%.0Lf  pieces where a shorter vector wins: %lld (max gain %.6Lf)\n",
           S3, floorl(S3), nonlongest, maxgainshort);
    return 0;
}
