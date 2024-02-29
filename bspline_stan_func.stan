// b-spline
    vector bspline(vector x, real xh, real delt) {
        
        int N = size(x);
        real xh1 = xh + 1 * delt;
        real xh2 = xh + 2 * delt;
        real xh3 = xh + 3 * delt;
        real xh4 = xh + 4 * delt;
        real u;
        vector[N] bs = rep_vector(0.0, N);
        
        for (i in 1:N) {
            if (xh <= x[i] && x[i] < xh1) {
                u = (x[i] - xh) / delt;
                bs[i] = 1.0 / 6 * u^3;
            }
            else if (xh1 <= x[i] && x[i] < xh2) {
                u = (x[i] - xh1) / delt;
                bs[i] = 1.0 / 6 * (1 + 3 * u + 3 * u^2 - 3 * u^3);
            }
            else if (xh2 <= x[i] && x[i] < xh3) {
                u = (x[i] - xh2) / delt;
                bs[i] = 1.0 / 6 * (4 - 6 * u^2 + 3 * u^3);
            }
            else if (xh3 <= x[i] && x[i] < xh4) {
                u = (x[i] - xh3) / delt;
                bs[i] = 1.0 / 6 * (1 - 3 * u + 3 * u^2 - u^3);
            }
            // else zeros [initialised as zeros]
        }   
        return bs;
    }
    
    // return Bspline design matrix, uniform grid
    matrix uniBspline(vector xx, int H) {
        // list of uniform B-splines over range [xx]
        // [COMPACT: extrapolation to zero]
        int N = size(xx);
        matrix[N, H] PHI;
        
        real delta;
        delta = (max(xx) - min(xx)) / (H + 3);
        
        for (h in 1:H) {
            PHI[,h] = bspline(xx, min(xx) + (h-1) * delta, delta);
        }
        return PHI;
    }

    // predict xp for Bsplines over predefined range [xx]
    matrix frac_uniBspline_pred(vector xp, vector xx, int H) {
        // i.e. the potential to extrapolate
        // as above [extrapolation to zero]
        int N = size(xp);
        matrix[N, H] PHI;
        real delta;
        delta = (max(xx) - min(xx)) / (H + 3);
        
        for (h in 1:H) {
            PHI[,h] = bspline(xp, min(xx) + (h-1) * delta, delta);
        }
        return PHI;
    }