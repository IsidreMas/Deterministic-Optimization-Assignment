# include <stdio.h>
# include <stdlib.h>

/*The number of values in x must match the dimension*/
float rosenbrock(float *x, unsigned d, float a, float b)
{
    float sum = 0, term1, term2;
    for(int i = 0; i<d-1; i++)
    {
        term1 = x[i+1]-x[i]*x[i];
        term2 = a - x[i];
        sum += b*term1*term1 + term2*term2;
    }
    return sum;
}

void gradient_rosenbrock(float *x, float *gradient, unsigned d, float a, float b)
{
    float term1, term2;
    for(int i = 0; i<d; i++)gradient[i]=0;
    for(int i = 0; i<d-1; i++)
    {
        term1 = 2*b*(x[i+1]-x[i]*x[i]);
        term2 = -2*term1*x[i]+2*(x[i]-a);
        if(i<d-1)
            gradient[i+1] += term1;
            gradient[i] += term2;
    }
}

float scalar_product(float *vec1, float *vec2, int dim)
{
    float sum = 0;
    for(int i = 0; i<dim; i++)sum+=vec1[i]*vec2[i];
    return sum;
}

void LM_d(float *x, float *d, float *gradient, float lambda) //2 dimensional -(H-lambda*I)^-1*gradient
{
    float A = (lambda+200)*(-lambda+400*x[1]-2)-400*(3*lambda+200)*x[0]*x[0];
    d[0] = -((-200-lambda)*gradient[0]-400*x[0]*gradient[1])/A;
    d[1] = -(-400*x[0]*gradient[0]+(-lambda-1200*x[0]*x[0]+400*x[1]-2)*gradient[1])/A;
}

void SteepestGradientDescent(float *seed, float *solution, unsigned dimension)
{
    FILE *file;
    file = fopen("../results/SteepestGradientDescent.dat","w");

    if(file == NULL)
    {
      printf("Error!");   
      exit(1);             
    }

    float *gradient,*xx, alphak, *d, fk, fk1, sigma = 1e-4, rho = 5e-1, min;
    unsigned iter = 0;
    if((gradient=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((d=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((xx=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    
    for(int i = 0; i < dimension; i++)solution[i]=seed[i];
    fprintf(file, "x_1\tx_2\tf(x_1,x_2)\n");
    do
    {
        min = rosenbrock(solution,2,1,100);
        fprintf(file, "%f\t%f\t%f\n", solution[0], solution[1], min);
        alphak = 1;
        gradient_rosenbrock(solution,gradient,dimension,1,100);
        fk =rosenbrock(solution,dimension,1,100);
        for(int i = 0; i<dimension; i++)
        {
            d[i] = -gradient[i];//define d
            xx[i] = solution[i];
            solution[i]+=d[i]*alphak;
        }
        fk1 = rosenbrock(solution,dimension,1,100);
        while(fk1>fk+sigma*alphak*scalar_product(gradient, d, dimension))
        {
            alphak*=rho;
            for(int i = 0; i<dimension; i++)solution[i]= xx[i] + d[i]*alphak;
            fk1 = rosenbrock(solution,dimension,1,100);
        }
        iter++;
    }while(rosenbrock(solution,2,1,100)<min);

    fclose(file);
    free(gradient);
    free(d);
    free(xx);
    printf("Iterations = %u\n", iter);
}

void ConjugateGradientMethod2(float *seed, float *solution, unsigned dimension)
{
    FILE *file;
    file = fopen("../results/ConjugateGradientMethod.dat","w");
    float betak = 0;

    if(file == NULL)
    {
      printf("Error!");   
      exit(1);             
    }

    float *gradient,*xx, alphak, *d,*dk,fk, fk1, sigma = 1e-1, rho = 5e-1, min;
    unsigned iter = 0;
    if((gradient=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((d=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((dk=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((xx=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    fprintf(file, "x_1\tx_2\tf(x_1,x_2)\n");
    for(int i = 0; i < dimension; i++)solution[i]=seed[i];
    do
    {
        min = rosenbrock(solution,2,1,100);
        fprintf(file, "%f\t%f\t%f\n", solution[0], solution[1], min);
        gradient_rosenbrock(solution,gradient,dimension,1,100);
        alphak=1;
        fk =rosenbrock(solution,dimension,1,100);
        for(int i = 0; i<dimension; i++)
        {
            d[i] = -gradient[i] + betak*d[i];//define d
            dk[i] = gradient[i];
            xx[i] = solution[i];
            solution[i]= xx[i] + d[i]*alphak;
        }
        fk1 = rosenbrock(solution,dimension,1,100);
        while(fk1>fk+sigma*alphak*scalar_product(gradient, d, dimension))
        {
            alphak*=rho;
            for(int i = 0; i<dimension; i++)solution[i]= xx[i] + d[i]*alphak;
            fk1 = rosenbrock(solution,dimension,1,100);
        }
        gradient_rosenbrock(solution,gradient,dimension,1,100);
        betak = scalar_product(gradient,gradient,dimension)/scalar_product(dk,dk,dimension);
        iter++;
    }while(rosenbrock(solution,2,1,100)<min);

    free(gradient);
    free(d);
    free(dk);
    free(xx);
    fclose(file);
    printf("Iterations = %u\n", iter);
}

void LevenvergMarquardt(float *seed, float *solution, unsigned dimension)
{
    FILE *file;
    file = fopen("../results/LevenbergMarquardt.dat","w");

    if(file == NULL)
    {
      printf("Error!");   
      exit(1);             
    }

    float *gradient,*d,*xx,lambda=1e-3,min;
    unsigned iter = 0;
    if((gradient=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((d=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    if((xx=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");

    for(int i = 0; i < dimension; i++)solution[i]=seed[i];
    fprintf(file, "x_1\tx_2\tf(x_1,x_2)\n");
    do
    {
        min = rosenbrock(solution,dimension,1,100);
        fprintf(file, "%f\t%f\t%f\n", solution[0], solution[1], min);
        gradient_rosenbrock(solution,gradient,dimension,1,100);
        LM_d(solution,d,gradient,lambda);
        for(int i = 0; i < dimension; i++)xx[i]=solution[i]+d[i];
        while(rosenbrock(xx,2,1,100)>min)
        {
            lambda*=5;
            LM_d(solution,d,gradient,lambda);
            for(int i = 0; i < dimension; i++)xx[i]=solution[i]+d[i];
        }
        lambda*=0.1;
        for(int i = 0; i < dimension; i++)solution[i]=xx[i];

        iter++;
    }while(rosenbrock(solution,dimension,1,100)<min);

    free(d);
    free(gradient);
    free(xx);
    fclose(file);
    printf("Iterations = %u\n", iter);


}



int main(int argc, char *argv[])
{
    float *x, *solution, a=1, b=100;
    unsigned dimension = 2;
    if((x=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for variables");
    if((solution=(float *)malloc(dimension*sizeof(float)))==NULL)
        printf("Couldn't allocate memory for gradient");
    
    x[0]=-1.5;
    x[1]=-1;
    if(argc>1)
        x[0] = atoi(argv[1]);
    if(argc>2)
        x[1] = atoi(argv[2]);

    gradient_rosenbrock(x,solution,dimension,1,100);
    printf("Gradient = (%f, %f)\n", solution[0], solution[1]);

    printf("f(-1.5,-1)=%f\n",rosenbrock(x,dimension,a,b));
    SteepestGradientDescent(x,solution,dimension);
    printf("seed(%f,%f)->f(%f,%f)=%f\n",x[0],x[1], solution[0], solution[1], rosenbrock(solution,dimension,a,b));
    ConjugateGradientMethod2(x,solution,dimension);
    printf("seed(%f,%f)->f(%f,%f)=%f\n",x[0],x[1], solution[0], solution[1], rosenbrock(solution,dimension,a,b));
    LevenvergMarquardt(x,solution,dimension);
    printf("seed(%f,%f)->f(%f,%f)=%f\n",x[0],x[1], solution[0], solution[1], rosenbrock(solution,dimension,a,b));
    return 0;
}