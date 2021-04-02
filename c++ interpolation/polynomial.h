#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include<vector>
#include<iostream>
#include<fstream>

class Polynomial
{
public:
    Polynomial();
    Polynomial(std::string& fileName);
    ~Polynomial();

    void save(std::string fileName);
    void load(std::string fileName);
    float eval(float x);

    friend
    std::ostream& operator<<(std::ostream& out,const Polynomial& x);
    friend
    std::istream& operator>>(std::istream& in,Polynomial& x);

private:
    std::vector<float> a;
    std::vector<float> xs;
    std::vector<float> ys;

    void coef();
    void parse(std::string stream,std::vector<float>& vector);
};

#endif // POLYNOMIAL_H
