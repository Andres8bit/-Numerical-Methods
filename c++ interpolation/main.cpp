#include <iostream>
#include"polynomial.h"
using namespace std;
bool isDone();
int main(int argc, char *argv[])
{

    string fileName = argv[1];//path of file containing points for intperpolation
    Polynomial p(fileName);
    string input = "";
    float x = 0;


    do{
        cout<<"x:";
        getline(cin,input);
        cout<<"input:"<<input<<endl;
        x = atof(input.c_str());
        cout<<endl<<"p("<<x<<")= "<<p.eval(x)<<endl;
    }while(!isDone());


    return 0;
}

bool isDone()
{
    string input = "";

    cout<<"quit?[q]:";
    getline(cin,input);

    if(input.at(0) == 'q' || input.at(0) == 'Q')
        return true;

    return false;
}
