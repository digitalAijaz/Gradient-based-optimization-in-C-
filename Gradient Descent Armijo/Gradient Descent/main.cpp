#include <iostream>
#include "Gradient.h"
#include <functional>

//Define the Object function
double ObjectFcn(std::vector<double>* funcLoc)
{
    //y = x^2
    double x = funcLoc->at(0);
    double y = funcLoc->at(1);
    //    double z = funcLoc ->at(2);
    //    return (x*x*x) + (2*(x*x)) - (2*x);
    //    return (x*x);
    //    return (x*x) + (2*(x*y)) + (y*y);
    //    return 4*x*x + y*y -2*x*y;
    //    return (x-7)*(x-7) + (y-2)*(y-2);
          return 100 * (y - x * x) * (y - x * x) + (1 - x) * (1 - x); 
    //    return 100*x*x*x*x + x*x -200*x*x*y -2*x +100*y*y +1;
    //    return (4*x*x) + x*(3*y) + (y*y) - (y*z) + (z*z) + (z);

}

using namespace std;

int main()
{
    //create a function pointer for our object function
    std::function<double(std::vector<double>*)> p_ObjectFcn{ ObjectFcn };

    //create a test instance of the qbGradient class
    Gradient solver;

    //Assign the object function
    solver.SetObjectFcn(p_ObjectFcn);

    //Set a start point
    std::vector<double> startPoint = {0.0,0.0 };
    solver.SetStartPoint(startPoint);

    //Set the maximum number of iterations
    solver.SetMaxIterations(5000);

    //set the step size
    solver.SetStepSize(0.001);

    //Call optimize
    std::vector<double> funcLoc;
    double funcVal;
    solver.Optimize(&funcLoc, &funcVal);

    //Output the result
   // std::cout<<"Function location: "<<funcLoc[0]<<std::endl;
    std::cout << "Function location: " << funcLoc[0] << " , " << funcLoc[1] << std::endl;
    //std::cout<<"Function location: "<<funcLoc[0]<<" , "<<funcLoc[1]<<" , "<<funcLoc[2]<<std::endl;
    std::cout << "Function value: " << funcVal << std::endl;

    return 0;
}
