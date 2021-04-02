#include "polynomial.h"

Polynomial::Polynomial()
{
    xs = ys = a = std::vector<float>();
}

Polynomial::Polynomial(std::string& fileName){
    load(fileName);
    coef();
}

Polynomial::~Polynomial(){
    xs = ys = a = std::vector<float>();
}

void Polynomial::load(std::string fileName){
    std::ifstream file;
    std::string tempX,tempY;
    file.open(fileName.c_str(),std::ios::binary);

    if(file.fail()){
        std::cout<<"failed to load"<<std::endl;
        file.close();
        return;
    }
    else{
        std::getline(file,tempX);
        std::getline(file,tempY);
        file.close();
    }
    parse(tempX,xs);
    parse(tempY,ys);
}

float Polynomial::eval(float x){
    unsigned int n = xs.size() - 1;
    float temp = 0;
        std::cout<<"call to eval :x = "<<x<<std::endl;
    for(int i = n; i >= 0 ;i--){
        temp = temp*(x - xs[i]) + a[i];
    }

    return temp;
}

//private:
void Polynomial::coef(){
    unsigned int n = ys.size();

    for(unsigned int i = 0; i < n;i++)
        a.push_back(ys[i]);

    for(unsigned int i = 1; i < n;i++)
        for(unsigned int j = n - 1; j >= i;j--)
             a[j] = (a[j] -a[j - 1])/(xs[j] - xs[j - i]);
}

void Polynomial::parse(std::string stream,std::vector<float>& vector){
    unsigned int i = 0;
    std::string val;

    while(stream.size()) {
        i = stream.find_first_of(" ");
        if(i == std::string::npos)
        {
            val = stream.substr(0,stream.size());
            stream.erase();

        }
        else
        {
            val = stream.substr(0,i);
            stream.erase(0,i + 1);
        }
        vector.push_back(std::atof(val.c_str()));
    }
}
