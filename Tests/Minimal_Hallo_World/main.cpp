#include <stdio.h>
#include <iostream>


int add_int(int a,int b)
{
	int sum = a+b;
	return sum;
}

int main()
{
    std::cout << "bonjour\n" << std::endl;
	printf("bonjour\n");
	
    // **************************
    // Addition de deux entiers
    // **************************
	int n=3,m=4;
	int p=n+m;
	printf("la somme de %d  et%d est égale à %d\n",n,m,p);
    std::cout << n << "  + "<< m<< " =  " << p << std::endl;

    // **************************
    // Avec une fonction 
    // **************************
	int q=add_int(n,m); 
	printf("la somme de %d  et%d est égale à %d\n",n,m,q);
	
	return 0;
}
