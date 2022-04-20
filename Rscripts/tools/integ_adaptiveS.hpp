#include <math.h>  // include file for fabs and sin
#include <stdio.h> // include file for printf


int main(float a, float b, float eps, float mrd){
  
  //float integ_adaptiveS(float a, float b, float eps, float mrd);        // compute integral of function(x)
  
    return adaptiveSimpsons(theta_spline, a, b, eps, mrd);
}

//
// Adaptive Simpson's Rule
//
float adaptiveSimpsons(float (*f)(float),   // ptr to function
                       float a, float b,  // interval [a,b]
                       float epsilon,  // error tolerance
                       int maxRecursionDepth) {   // recursion cap        
  float c = (a + b)/2, h = b - a;                                                                  
  float fa = f(a), fb = f(b), fc = f(c);                                                           
  float S = (h/6)*(fa + 4*fc + fb);                                                                
  return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}                                                                                                   

//
// Recursive auxiliary function for adaptiveSimpsons() function below
//                                                                                                 
float adaptiveSimpsonsAux(float (*f)(float), float a, float b, float epsilon,                 
                          float S, float fa, float fb, float fc, int bottom) {                 
  float c = (a + b)/2, h = b - a;                                                                  
  float d = (a + c)/2, e = (c + b)/2;                                                              
  float fd = f(d), fe = f(e);                                                                      
  float Sleft = (h/12)*(fa + 4*fd + fc);                                                           
  float Sright = (h/12)*(fc + 4*fe + fb);                                                          
  float S2 = Sleft + Sright;                                                                       
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)   // magic 15 comes from error analysis                                       
    return S2 + (S2 - S)/15;                                                                        
  return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
    adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}

float theta_spline(float t, float theta0, float nu) {
  
  float theta = theta0 * exp(-nu * t)
}