#include <math_ext.h>

///////////////////////////////////////////////////////////////////////////////
//     Math/Stat functions                                                   //
///////////////////////////////////////////////////////////////////////////////

// Returns the vector of differences between first and
// last elements of the windows of size n in x
// from igraph_vector_t to  std::vector
int rolling_difference_igraph(igraph_vector_t &x, std::vector<double> &out, int n){

    int len_out = igraph_vector_size(&x) - n + 1;
    out.resize(len_out);

    // for each window start
    for(int ind=0; ind < len_out; ind++) {
        // difference is between first and last elements
        // of the window
        out[ind] = VECTOR(x)[ind+n-1] - VECTOR(x)[ind];
    }

    return 0;
}


// Returns the vector of differences between first and
// last elements of the windows of size n in x
// from std::vector to  std::vector
// Nope nope nope TODO
int rolling_difference(std::vector<double> &x, std::vector<double> &out, int n){

    int len_out = x.size() - n;
    out.resize(len_out);

    // for each window start
    for(int ind=0; ind < len_out; ind++) {
        // difference is between first and last elements
        // of the window
        out[ind] = x[ind + n] - x[ind];
    }

    return 0;
}


// Median of a vector
double median(std::vector<double> v){

    int size = v.size();
    if(size == 0){
        return NAN;
    }

    else if(size == 1){
        return v[0];
    }

    else{
        if(size % 2 == 0){
            // even length vector, average of middle 2 numbers
            double o0;
            double o1;
            std::nth_element(v.begin(), v.begin() + size/2, v.end());
            o1 = v[size/2];
            o0 = *std::max_element(v.begin(), v.begin() + size/2 -1);
            return (o0 + o1)/2;
        }
        else{
            // odd length, middle number
            std::nth_element(v.begin(), v.begin() + (size-1)/2, v.end());
            return v[(size-1)/2];
        }
    }
}

// Mean/average of a vector
double mean(std::vector<double> v){
    double mean = 0.0;
    int n = v.size();
    for(int i=0; i<n; i++){
        mean = mean + v[i];
    }
    return mean/n;
}

// Standard deviation of a vector
double stddev(std::vector<double> v, double dof){
    double n = v.size();
    double v_bar = mean(v);
    double var = 0;

    for(int i=0; i<n; i++){
        var = var + pow( (v[i] - v_bar), 2.0);
    }

    var = var / (n-dof);
    return sqrt(var);
}

// get the exponent to pow value to make a double an int
// Stephen Grady
int get_precision(double k){

    int exponent=0;
    bool doneWith0=false;

    std::ostringstream convert;
    convert<<k;
    std::string str=convert.str();

    for (int i = str.size()-1; i > -1; --i){
        if(str[i]=='.'){
            return exponent;
        }
        // first nonzero value, start counting places
        if(str[i]!='0'){
            doneWith0=true;
        }
        if(doneWith0){
            exponent++;
        }
    }
    return 0;
}

// Range from l to u, incrementing by increment
std::vector<double> range(double l, double u, double increment){
    // doubleing point arithmetic
    double precision = pow(10, get_precision(increment));

    int int_increment = static_cast<int>(increment*precision);
    int int_l = static_cast<int>(l*precision);
    int int_u = static_cast<int>(u*precision);

    std::vector<double> out;

    for(int int_t=int_l; int_t<=int_u; int_t+=int_increment){
        out.push_back(int_t/precision);
    }

    return out;
}

// Empirical Cumulative Distribution Function (ecdf)
std::vector<double> ecdf(std::vector<double> x, std::vector<double> t){
    // Evaluate the ecdf of x at points t
    // Sorted x, so for each x_i in x, F(x_i) = index of x_i +1 / Number observations
    // Instead we are evaluating at t, F(t_i) = (number of x_i smaller than t_i) / Number observations

    if( !std::is_sorted(t.begin(), t.end()) ){
        throw std::invalid_argument("ecdf evaluate at points is not sorted");
    }

    if( !std::is_sorted(x.begin(), x.end()) ){
        std::sort(x.begin(), x.end());
    }

    int n = x.size();   // Number observations
    int N = t.size();   // Number places to evaluate F

    std::vector<double> F(N);

    // j = current index in x
    // i = current index in t
    // k = placeholder index of x

    int k=0, j=0, i = 0;

    for(i=0; i<N; i++){

        if(j<n){
            for(j=k; j<n; j++){
                if(x[j] > t[i]){
                    break;
                }
                else{
                    k=j;
                }
            }
            F[i] = ( (double) j ) / (double) n;
        }
        else{
            F[i] = 1;
        }

        //std::cout << "ecdf\ti= " << i << "\tj= " << j << "\tk= " << k << "\t\tt[i]= " << t[i] << "\t\tx[j]= " << x[j] << "\t\tF[i]= " << F[i] << std::endl;
    }
    return F;
}

double sign(double x){
    if(x == 0){
        return 1.0;
    }
    else{
        return std::copysign(1.0, x);
    }
}

// These two functions return area under distribuition between x1 and x2, x1 < x2
double poisson(double x1, double x2){
    // p(s) ds = exp(-s) ds
    return exp(-x1) - exp(-x2);
}

double goe(double x1, double x2){
    // p(s) d(s) = 0.5 * M_PI * x * exp( -1.0 * M_PI * pow(x, 2.0) / 4.0) ds
    return exp( -1.0 * M_PI * pow(x1, 2.0) / 4.0) - exp( -1.0 * M_PI * pow(x2, 2.0) / 4.0);
}

double fisher_transform(double r, int n){
    return atanh(r) + r / (2.0 * (double(n) - 1.0));
}