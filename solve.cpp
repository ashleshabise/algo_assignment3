#include<stdio.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#define a 16
#define SIZE 100 
#define TRUE 1
#define FALSE 0

 fGamma(int K, float q, char *partitionOutputFile, char *inputFileGamma, char *fileForR) {
 
 char *str;
 int  gammaArray[K], i=0;
 // Read gammas from inputFileGamma
 FILE gammaFile = fopen(inputFileGamma,"r");
	while (fgets(str,SIZE,gammaFile) != NULL) {
		sscanf(str,"%f",&gamma);
		gammaArray[i++] = gamma;
	}

	// Read file for R
	FILE rFile = fopen(fileForR, "r");

// 
 float rArray[3];
	while (fgets(str,SIZE,rFile) != NULL) {
		sscanf(str,"%f %f %f",&x1,&x2,&x3);
		rArray[0] = x1;
		rArray[1] = x2;
		rArray[2] = x3;
		count++;
	}



	
	FILE pOutputFile = fopen(partitionOutputFile, "r");
	while (fgets(str,SIZE,pOutputFile) != NULL) {
		sscanf(str,"%f %f %f",&x1,&x2,&x2);

	}
 }



int summation(int l[] [], int y[],int start, int end, int i) {
	int sum = 0;
	for(int j=start; j<=end; j++) {
		sum+=l[i][j]*y[j]
	}
	return sum;	
}

int* LUP-SOLVE(int L[][], int U[][],int Pi[],int b[], int rowsL) {
	

	int n=rowsL;
	int P[n][n], y[n], x[n];

	// Initialize P matrix
	for ( int i=0; i < n; i++) {
		for ( int j=0; j < n; j++) {
			if ( j == Pi[i])
				P[i][j] = 1;
			else
				P[i][j] = 0;
		}
	}
	// Calculate Pb from P and b
	for ( int i=0; i < n; i++) {
		for ( int j=0; j < n; j++) {
			Pb[i] += P[i][j]*b[j];
		}
	}

	// Initialize y and x
	for (int i=0; i < n; i++ )
		y[i] = 0,x[i] = 0;

	for(int i=1; i<n; i++)
		y[i] = Pb[i]- summation(L,y,1,i-1);

	for(int i=n;i>=0;i--) 
		x[i]=(y[i]-summation(U,x,i+1,n))/U[i][i];
	
	int *result = x;
	// returned a pointer to the array x.
	return result;	
}
LUP-Decomposition(int A[][]) {
	int n=rowsA, temp;
	int Pi[n];

	for(int i=1;i<=n;i++)
		Pi[i]=i;

	int p, k1;
	for(int k=1; k<n; k++) {
		p=0;
		for(int i=k; i<n; i++) {
			if(abs(A[i][k]) > p) {
				p=abs(A[i][k]);
				k1=i;
			}
		}
		if (p==0)
			cout<<"\n Singular matrix.";

		// Swap 
		temp=Pi[k];
		Pi[k]=Pi[k1];
		Pi[k1]=temp;

		for(int i=1; i<n; i++) {
			temp=A[k][i];
			A[k][i]=A[k1][i];
			A[k1][i]=temp;
		}

		for(int i=k+1; i<=n; i++) {
			A[i][k]/= A[k][k];
			for(int i=1; i<n; i++)
				A[i][j]-= (A[i][k] * A[k][j]);
		}
		
	}
}

// Returns a pointer to solution array.
int* solve-ls(int A[][], int b[], int ) {
	LUP-Decomposition();
	int* result = LUP-SOLVE();
	return result;
}

/* Code for SOLVE function goes here */
bool solve(int L, char* inputFileR, char* inputFileP, char* inputFileQ, char* inputFileGamma, char* outputFile) {
	
}


/*--------------------------------------------------------------------------------------------*/

/* Code for partitionAll + utility functions starts here */

// A linked list where each node represents a point r[]
struct node {
	float r[3];
	struct node* next;
};
void push(struct node** head, float r[3]) {
	struct node* newNode = (struct node*)(malloc(sizeof(struct node)));
	newNode->next = NULL;
	int i;
	// copy r[i] into newNode->r[i];
	for (i=0; i<3; i++) newNode->r[i] = r[i];
	if ( *head == NULL) {
		*head = newNode;
		return;
	}
	newNode->next = *head;
	*head = newNode;
}

bool partition(float r[3],int L, int lMatrix[L][L][L], struct node* pointList[L][L][L]) {
	//validity check for r[3]
	fflush(stdout);
	printf("Entering partition\n");
	int i;
	for(i=0;i<3;i++) { 
		if(!(r[i]>=0 && r[i] <= a)) {  
			printf(" Input Invalid");
					return FALSE;
		}
		printf("%f ", r[i]);
	}
	printf("\n");
	//validity check for L
	if(!(L>0))
		printf(" Input Invalid");

	FILE* fPartition;
	fPartition=fopen("partitionFile", "a");
	if(fPartition==NULL) {
		printf("Couldn't open partition file");
		return FALSE;
	}
	float x[3];
	printf("opened partition file\n");
	int l[3],binary[3], bString=0, temp;
	for (i=0; i<3; i++)	{
	printf("generating binary string: ");
	 x[i]=(r[i]*L)/a;
	 l[i] = (int)x[i];
	 // for x1, temp=1. binary = 0 | 1 = 1
	 // for x2, temp=0, binary =  1<<1 = 10 | 0 = 10
	 // for x3, temp=1, binary  = 10<<1  = 100 | 1 = 101
	 printf("%f %d", x[i], l[i]);
	 binary[i] =  (x[i]-l[i]) ? 0 : 1;
	 printf("binary[%d] : %d \n",i, binary[i]);
	 bString = bString<<1 | binary[i];
	} 
	printf("\n");
	// write 000 always. 
	fprintf(fPartition,"%d %d %d %f %f %f\n", l[0], l[1], l[2], r[0],r[1],r[2]);
	// update the lMatrix as well. 
	lMatrix[l[0]][l[1]][l[2]]++;
	int b1,b2,b3;
	int done[8]={0,0,0,0,0,0,0,0},result;
	for ( i=1; i <= bString; i++)   {
	printf("writing l[] into partition file\n");
		result = bString & i;
		if (result ) {
			b1=result&4 ? 1:0;
			b2=result&2 ? 1:0;
			b3=result&1 ? 1:0;
			// update the lMatrix
			if (!done[result]) {
				fprintf(fPartition,"%d %d %d %f %f %f\n", l[0]-b1, l[1]-b2, l[2]-b3, r[0],r[1],r[2]);
				lMatrix[l[0]-b1][l[1]-b2][l[2]-b3]++;
				// push the point corresponding to l[] to the pointList 
				push(&pointList[l[0]-b1][l[1]-b2][l[2]-b3],r);
				done[result]=1;
			}
			printf("%d\n", done[result]);
		}
	}
	fclose(fPartition);
	return TRUE;
}
bool partitionAll(int L, char* inputFileName, char* outputFileName) {
	printf("Entering partitionall\n");
	float r[3];
	float **R;
	// Declaring helper variables 
	int i,j,k;
	int a1,a2,a3;
	float x1,x2,x3, temp1, temp2, temp3;
	int lMatrix[L][L][L],K=0;
	// initialise every element of lMatrix to 0. 
	for ( i=0; i<L; i++) 
		for ( j=0; j<L; j++) 
			for ( k=0; k<L; k++) 
				lMatrix[i][j][k]=0;

	FILE *fin, *fout;
	fin = fopen(inputFileName,"r");
	if (fin != NULL) {
		char str[SIZE];
		while (fgets(str,SIZE,fin) != NULL) {
			printf("getting data from input file\n");
			/*
			temp1 = x1;
			temp2 = x2;
			temp3 = x3;
			*/
			sscanf(str, "%f %f %f", &x1, &x2, &x3);
			printf("%f %f %f\n", x1,x2,x3);
			if (x1>a || x1<0 /*|| x1==temp1*/ ) { 
				printf("Input Failed\n");
				return FALSE;
			}
			if (x2>a || x2<0 /*|| x2==temp2*/) {
				printf("Input Failed\n");
				return FALSE;
			}
			if (x3>a || x3<0 /*|| x3==temp3*/) {
				printf("Input Failed\n");
				return FALSE;
			}
			r[0] = x1;
			r[1] = x2;
			r[2] = x3;
			partition(r,L,lMatrix);
			// so we have K points. 
			K++;
		}
	}
	else {
		printf("Could not open the input File  \n");
		return 1;
	}
	fclose(fin);
	fout = fopen(outputFileName, "w");
	FILE *fPartition ;
	char str[SIZE];
	for ( i=0; i<L; i++) {
		for ( j=0; j<L; j++) {
			for ( k=0; k<L; k++) {
				fprintf(fout,"%d %d %d\n", i, j, k);
				fprintf(fout,"%d\n",lMatrix[i][j][k]);
				// open the partition file
				fPartition = fopen("partitionFile", "r");
				while (fgets(str,SIZE,fPartition) != NULL) {
					sscanf(str,"%d %d %d %f %f %f",&a1,&a2,&a3,&x1,&x2,&x2);
					if (i==a1 && j==a2 && k==a3)
						fprintf(fout,"%f %f %f\n",x1,x2,x3);
				}
				// close the partition file.
				fclose(fPartition);
			}
		}
	}
	fclose(fout);
	return TRUE;
}

/*----------------------------------------------------------------------------------*/

/* Code for BUILD-MATRIX function starts here */

bool MATRIX(int L, int l[3], char* fileForP, char* fileForQ, char* fileForMatrix) {
	float x[3], temp1, temp2, temp3;
	float x1,x2,x3;
	FILE *fp, *fq, *fMatrix;
	int N=0, M=0,i,j;
	float P[100][3],Q[100][3];
	// Initialise P and Q to 0
	for (i=0; i < 100; i++) {
		for ( j=0;j<3; j++) {
			P[i][j] = Q[i][j] = 0;
		}
	}
	fp = fopen(fileForP,"r");
	if (fp != NULL) {
		char str[SIZE];
		// Checking validity of input P and pushing valid data into P matrix
		while (fgets(str,SIZE,fp) != NULL) {
			sscanf(str, "%f%f%f", &x[0], &x[1], &x[2]);
			int flag = 0;
			for ( i=0; i<3; i++) {
			// check for 0<= x[i] <= a
				if (x[i]>a || x[i]<0) {
					printf("initial condition : Invalid Input for P");
					return FALSE;
				}
				if (!( (((a*l[i])/L == x[i]) ||(x[i] == (a*(l[i]+1)/L)))   &&   ( ((a*l[(i+1)%3])/L <= x[(i+1)%3]) && (x[(i+1)%3]<= (a*(l[(i+1)%3]+1)/L)) )   &&   ( ((a*l[(i+2)%3])/L <= x[(i+2)%3]) &&  (x[(i+2)%3]<= (a*(l[(i+2)%3]+1))/L) ))) {
					flag++;
				}
			}
			/*
				If the condition fails for all three co-ordinates,
				the point is not a boundary point
			*/
			if (flag==3) {
				printf("Invalid Input for P\n");
				return FALSE;
			}
			P[N][0] = x[0];
			P[N][1] = x[1];
			P[N][2] = x[2];
 			/* This will count no. of inputs in file for P*/
			N++;
		}
	} else {
		printf("Could not open the input File  \n");
		return 1;
	}
	fclose(fp);
	fq = fopen(fileForQ,"r");
	// Checking validity of input Q and pushing valid data into Q Matrix
	if (fq != NULL) {
		char str2[SIZE];
		while (fgets(str2,SIZE,fq) != NULL) {
			sscanf(str2, "%f%f%f", &x[0], &x[1], &x[2]);
			// Calculate the radius of the sphere. 
			float sum= powf(2*a/(float)L,2);
			// Check if the point lies on the sphere.
			for(i=0; i<3 ; i++)
				sum-= powf(x[i]-((2*l[i]+1)*a)/(2*L), 2);

			// Decide the validity accordingly
			if(sum) {
					printf("Input data for Q is Invalid");
					return false; 
				}
			Q[M][0] = x[0];
			Q[M][1] = x[1];
			Q[M][2] = x[2];
 			/* This will count no. of inputs in file for Q*/
			M++;
		}
 } else {
		printf("Could not open the input File  \n");
		return 1;
	}
	fclose(fq);
	// Construction of matrix A 
	float A[M][N]; 
	for ( i =0; i < M; i++) {
		for ( j =0; j < N; j++) {
			float temp;
			// sqrt( x^2 + y^2 + z^2 ) 
			temp = sqrt( powf((Q[i][0] - P[j][0]), 2) + powf((Q[i][1] - P[j][1]), 2) + powf((Q[i][2] - P[j][2]), 2));
			// here ?? 
			// so we need to multiply each element of a[i][j] with gamma  and add .
			// umm.. kaise ?? i mean.  
			A[i][j]= temp ? (sinf(temp))/temp : 1;
		}
	}
	// Finding the Product Matrix of A and Transpose of A 
	float result[N][N];
		for ( i =0; i < N; i++) {
			for ( j =0; j <= i; j++) {
				result[i][j] = 0;
				for ( int k=0; k<N; k++) {
					result[i][j]+= A[k][i] * A[k][j];
				}
				result[j][i] = result[i][j];
			}
		}
	fMatrix = fopen(fileForMatrix, "w");
		for ( i =0; i < M; i++) {
			for ( j =0; j < M; j++) {
					fprintf(fMatrix, "%d %d %f\n", i, j, result[i][j]);
				}
		}
	fclose(fMatrix);
}
				
/*  Main for build-matrix program

int main( int argc, char* argv[]) {
// check that whether file names are provided or not at the command line.
	if(argc!=4) { 
		printf("File Names are not provided at command line");
		return 1;
	}
	int L;
	do { 
		printf("Enter a Positive Integer: ");
		scanf("%d",&L);
	}
	while(L<=0);
	int l[3];
	do { 
		printf("\nEnter 3 integers less than %d: ",L);
		scanf("%d%d%d",&l[0],&l[1],&l[2]);
	}
	while(!(l[0]<L && l[1]<L && l[2]<L && l[0]>=0 && l[1]>=0 && l[2]>=0));
	bool result=MATRIX(L,l,argv[1],argv[2],argv[3]);
	switch(result) {
		case true: printf("\nAlgorithm is implemented Successfully\n");
								break;
		case false: printf("\nAlgorithm is implemented Successfully\n");
								break;
	}
	return 0;
}
*/


