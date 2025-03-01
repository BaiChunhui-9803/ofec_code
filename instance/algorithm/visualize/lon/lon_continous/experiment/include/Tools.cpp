#include "Tools.h"
#include <sstream>


using namespace std;
long long hashing_value(double x, int digits) {
    double a;
    double prec;
    long long b;

    if (digits == 0) {
        a = ceil(x);
    }
    else {
        prec = pow(10.0, -digits);
        a = x / prec;
        a = round(a);
    }
    long long from = -2147483647;
    long long to = 2147483647;


    if (a > to) {
        cout << " warning (hashing_value): out_of_range\t" << a<< endl;
        return to;
    }
    if (a < from) {
        cout << " warning (hasshin_value): out_of_range\t" << a<< endl;
        return from;
    }
    b = (long long)a;

    return b;
}

std::string toHashStr(const std::vector<double>& x, int digit) {

    std::stringstream out;
    for (int i = 0; i < x.size(); i++) {
        out << std::setfill('0') << setw(8) << hex << hashing_value(x[i],digit);
    }
    return std::move(out.str());
}

std::string toHashStrFit(double fit, int digit) {
    return to_string(hashing_value(fit, digit));

    std::stringstream out;
    out << setw(11) << dec << hashing_value(fit, digit);
    return std::move(out.str());
}

bool loadConfigurationFile(const string name, vector<string>& func,
    vector<int>& nvar, vector<int>& bounded,
    vector<double>& best,
    string& step_mode, vector<double>& beta,
    vector<double>& factor, int& opt_digits,
    int& hash_digits, int& iter, int& runs,
    string& init_mode, string& output_dir) {
    fstream fin;
    int nfunc, nfac;
    string cmd;

    fin.open(name, ios::in);
    if (!fin.is_open()) {
        cout << name << " not found !!!\n";
        return false;
    }

    func.clear();
    nvar.clear();
    best.clear();
    beta.clear();
    factor.clear();

    fin >> cmd >> nfunc >> cmd;
    func.resize(nfunc);
    for (int i = 0; i < nfunc; i++)
        fin >> func[i];

    nvar.resize(nfunc);
    fin >> cmd;
    for (int i = 0; i < nfunc; i++)
        fin >> nvar[i];

    bounded.resize(nfunc);
    fin >> cmd;
    for (int i = 0; i < nfunc; i++)
        fin >> bounded[i];

    best.resize(nfunc);
    fin >> cmd;
    for (int i = 0; i < nfunc; i++)
        fin >> best[i];

    fin >> cmd >> step_mode >> cmd;
    beta.resize(nfunc);
    for (int i = 0; i < nfunc; i++)
        fin >> beta[i];

    fin >> cmd >> nfac >> cmd;
    factor.resize(nfac);
    for (int i = 0; i < nfac; i++)
        fin >> factor[i];

    fin >> cmd >> opt_digits;
    fin >> cmd >> hash_digits;
    fin >> cmd >> iter;
    fin >> cmd >> runs;
    fin >> cmd >> init_mode;
    fin >> cmd >> output_dir;

    fin.close();
    cout << " " << name << " was loaded!!!\n";

    return true;
}

bool displayConfiguration(const vector<string>& func,
    const vector<int>& nvar, const vector<int>& bounded,
    const vector<double>& best,
    const string& step_mode, const vector<double>& beta,
    const vector<double>& factor, const int& opt_digits,
    const int& hash_digits, const int& iter, const int& runs,
    const string& init_mode, const string& output_dir) {
    cout << "functions: ";
    for (unsigned int i = 0; i < func.size(); i++)
        cout << func[i] << " ";
    cout << "\n";

    cout << "nvar: ";
    for (unsigned int i = 0; i < nvar.size(); i++)
        cout << nvar[i] << " ";
    cout << "\n";

    cout << "bounded: ";
    for (unsigned int i = 0; i < bounded.size(); i++)
        cout << bounded[i] << " ";
    cout << "\n";

    cout << "best: ";
    for (unsigned int i = 0; i < best.size(); i++)
        cout << best[i] << " ";
    cout << "\n";

    cout << "beta: ";
    for (unsigned int i = 0; i < beta.size(); i++)
        cout << beta[i] << " ";
    cout << "\n";

    cout << "factor: ";
    for (unsigned int i = 0; i < factor.size(); i++)
        cout << factor[i] << " ";
    cout << "\n";

    cout << "opt_digits: " << opt_digits << "\n";
    cout << "hash_digits: " << hash_digits << "\n";
    cout << "iter: " << iter << "\n";
    cout << "runs: " << runs << "\n";
    cout << "init_mode: " << init_mode << "\n";
    cout << "output_dir: " << output_dir << "\n";

    return true;
}