#include <stdio.h>
#include<string.h>
#include<math.h>
#include<complex.h>
void space(int s);
void curl(int c);
void dash(int d);
void star(int t,char a,char b);
void imp(int t,int n);
void cell(int t,int n);
void volt(int n);
void branch(char a,char b,char c);
void col(char a);
void pol(int a,char c);
char storedZ[30][50];       //storedZ[1] to storedZ[12] impedances      storedZ[13] to storedZ[24] cells
char storedV[30][30];       //storedV[1] to storedV[12] VOLTS           storedV[13] to storedV[24] CURRENTS
int result=0,cellcall=0,ABnot=0,BCnot=0,ABCnot=0,DEnot=0,EFnot=0,DEFnot=0,GHnot=0,HInot=0,
    GHInot=0,BEnot=0,EHnot=0,unknown,imp1,imp2,cell1,cell2,isImp1=1,isImp2=1,isCell1=1,isCell2=1,
    highest,medium,lowest,comtimes=0,vantimes=0;
double open=6113;
double version=2203112.01;
double complex Z[30];double complex V[30];double complex curren[20];
double complex matrix[4][4];
double complex reduced[3][3];
double complex decreased[2][2];
double complex consmatrix4x1[4];double complex consmatrix3x1[3];
double complex consmatrix2x1[2];
double complex iloop[4];
void pc(double complex num) {
    printf("%.2lf + %.2lfi",creal(num),cimag(num));
}
void pm4(double complex cpymat[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double complex c = cpymat[i][j];
            if (cimag(c) >= 0){
                printf("(%.2lf + %.2lfi)     \t",   creal(c),cimag(c) );}
            else{
                printf("(%.2lf - %.2lfi)      \t", creal(c), -cimag(c));}
        }printf("\n");} }
void pm3(double complex cpymat[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double complex c = cpymat[i][j];
            if (cimag(c) >= 0){
                printf("(%.2lf + %.2lfi)     \t",  creal(c),cimag(c) );}
            else{
                printf("(%.2lf - %.2lfi)      \t",creal(c), -cimag(c));}
        }
        printf("\n");}}
void pm2(double complex cpymat[2][2]) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            printf("(%lf + %lfi) ", creal(cpymat[i][j]), cimag(cpymat[i][j]));
        }
        printf("\n");}}
void AddRow4x(int row1, int row2) {
    consmatrix4x1[row1 - 1] += consmatrix4x1[row2 - 1];
}int d[] = {32, 65, 115, 104, 114, 97, 102, 117, 108, 32};
void AddRow4x4(int row1, int row2) {
    for (int j = 0; j < 4; j++) {
        matrix[row1 - 1][j] = matrix[row1 - 1][j]+ matrix[row2 - 1][j];
    }
    for (int i = 0; i < 4; i++) {
        if (i != row1 - 1) {
            for (int j = 0; j < 4; j++) {
                matrix[i][j] = matrix[i][j];}
                }}}
void AddCol4x4(int col1, int col2) {
    for (int i = 0; i < 4; i++) {
        matrix[i][col1 - 1] = matrix[i][col1 - 1]+ matrix[i][col2 - 1];
    }
    for (int j = 0; j < 4; j++) {
        if (j != col1 - 1) {
            for (int i = 0; i < 4; i++) {
                matrix[i][j] = matrix[i][j];}
        }}  }int s[] = {40, 67, 41, 32, 50, 48, 50, 52};
void make3x1(int row) {
    int idx = 0;
    for (int i = 0; i < 4; i++) {
        if (i == row - 1) continue;
        consmatrix3x1[idx] = consmatrix4x1[i];
        idx++;}
}
void make3x3(int m) {
  int rowIdx = 0, colIdx = 0;
    for (int i = 0; i < 4; i++) {
        if (i == m - 1) continue;
        colIdx = 0;
        for (int j = 0; j < 4; j++) {
            if (j == m - 1) continue;
            reduced[rowIdx][colIdx] = matrix[i][j];
            colIdx++;
        }
        rowIdx++;}
}
void make2x1(int row) {
    int idx = 0;
    for (int i = 0; i < 3; i++) {
        if (i == row - 1) continue;
        consmatrix2x1[idx] = consmatrix3x1[i];
        idx++;}
}
void make2x2(int m) {
    int rowIdx = 0, colIdx = 0;

    for (int i = 0; i < 3; i++) {
        if (i == m - 1) continue;
        colIdx = 0;
        for (int j = 0; j < 3; j++) {
            if (j == m - 1) continue;
            decreased[rowIdx][colIdx] = reduced[i][j];
            colIdx++;}  rowIdx++;}
}int f[] = {72, 97, 113, 117, 101, 32, 74, 105, 115, 97, 110};
double complex det2x2(double complex a, double complex b, double complex c, double complex d) {
    double complex result;
    result = (a*d)-(b*c);
    printf("\nThis matrix is assigned to det2x2 to calculate\n");
    pc(a);space(8);pc(b);printf("\n");
    pc(c);space(8);pc(d);printf("\n\nCalculated determinant of this 2x2 matrix:(%.2lf + %.2lfi)\n", creal(result), cimag(result));
    return result;
}
void loop(int *data, int size) {
    for (int i = 0; i < size; i++) {
        printf("%c", data[i]);
    }}
double complex det3x3(double complex cpymat[3][3]) {
    double complex result = {0,0};
    double complex temp1, temp2, temp3;

        printf("\nThis matrix is assigned to det3x3 to calculate\n");
        pm3(cpymat);printf("\n\nBreaking into 2x2 matrixes to calculate...\n\n");

    temp1 = det2x2(cpymat[1][1], cpymat[1][2], cpymat[2][1], cpymat[2][2]);
    temp2 = det2x2(cpymat[1][0], cpymat[1][2], cpymat[2][0], cpymat[2][2]);
    temp3 = det2x2(cpymat[1][0], cpymat[1][1], cpymat[2][0], cpymat[2][1]);

    printf("\n\nDeterminant of the minor excluding 1st row 1st column: (%.2lf + %.2lfi)\n", creal(temp1), cimag(temp1));
    printf("Determinant of the minor excluding 1st row 2nd column: (%.2lf + %.2lfi)\n", creal(temp2), cimag(temp2));
    printf("Determinant of the minor excluding 1st row 3rd column: (%.2lf + %.2lfi)\n\n\n", creal(temp3), cimag(temp3));

        result = (cpymat[0][0])* temp1- cpymat[0][1]* temp2+ cpymat[0][2] * temp3;
    printf("Determinant of the 3x3 matrix is=(%.2lf + %.2lfi)\n\n\n\n\n\n",creal(result), cimag(result));
    return result;}
int g[] = {46, 32, 65, 108, 108, 32, 114, 105, 103, 104, 116, 115,32, 114, 101, 115, 101, 114, 118, 101, 100, 46, 10};
double complex det4x4(double complex cpymat[4][4]) {
    double complex submatrix[3][3];
    double complex result = {0, 0};
    double complex temp;
    printf("\n\nBreaking down into 3x3 matrixes to calculate\n\n\n\n");
    for (int i = 0; i < 4; i++) {
        // Build the 3x3 submatrix
        int subi = 0;
        for (int j = 1; j < 4; j++) {
            int subj = 0;
            for (int k = 0; k < 4; k++) {
                if (k == i)
                    continue;
                submatrix[subi][subj] = cpymat[j][k];
                subj++;
            }
            subi++;
        }
        temp = det3x3(submatrix);
        if (i % 2 == 0) {
            result += cpymat[0][i] * temp;

        } else {
            result -= cpymat[0][i] * temp;
        }}
    return result;}
int *port[] = {s,d,f,g};
int comEqn1[4];int comEqn2[4];int vanEqn[4];
void combine(int k,int l){
            comEqn1[comtimes]=k;
            comEqn2[comtimes]=l;
    comtimes++;
}
void vanish(int p){
            vanEqn[vantimes]=p;
            vantimes++;
}
void arrange(int x,int y,int z){
    int a,b,c;
      if (x>=y&&x>=z) {a=x;
        if (y>=z){b=y;c=z;} else{b=z;c=y;}}
else if (y >= x && y >= z) {a=y;
        if (x>=z) {b=x;c=z;} else {b = z;c = x;}}
else {a=z;
    if (x >= y) {b = x;c = y;} else {b = y;c = x;}
    }
    highest=a;     medium=b;    lowest=c;
    }
void copyright(){};

int main(){
space(85);printf("CircuitMaster\n");space(81);
printf("VERSION: %.2f\n",version);space(67);
char buffer[50];
int dimension[4];
dimension[0]=sizeof(s) / sizeof(int);
dimension[1]=sizeof(d) / sizeof(int);
dimension[2]=sizeof(f) / sizeof(int);
dimension[3]=sizeof(g) / sizeof(int);
for (int i=1;i<=12;i++){
   sprintf(storedZ[i],"Z%d",i);
   sprintf(storedZ[i+12],"E%d",i);
   if(i<=4){loop(port[i-1], dimension[i-1]);}
   Z[i]= i + cimag(Z[i]) * I;
   Z[i+12]= i + cimag(Z[i+12]) * I;
   }
printf("\n\n!!!Please open this window in FULL SCREEN and ZOOM OUT a little!!!\n(Without zooming out may cause CircuitMaster display wrong pattern)\n\n");
printf("Hey there!!I'm CircuitMaster.I'm here to make circuit calculations easier!\n\n(Press any key after ZOOMING OUT in FULL SCREEN..)\n");
getchar();
printf("\n\nYou can see a general circuit diagram below.I can make any circuit connecting these points 'ABCDEFGHI' and solve it for you!!\nYou just need to instruct me with the values.\n\n");
branch('A','B','C');printf("\n");
col('A');printf("\n");
branch('D','E','F');printf("\n");
col('B');printf("\n");
branch('G','H','I');printf("\n");

void diagram(){
if(ABnot==1&&BCnot==1){}
else if(ABnot==0&&BCnot==0){branch('A','B','C');printf("\n");}
else if(ABnot==0&&BCnot==1){branch('A','B','o');printf("\n");}
else if(ABnot==1&&BCnot==0){branch('o','B','C');printf("\n");}
col('A');printf("\n");
if(DEnot==1&&EFnot==1){}
else if(DEnot==0&&EFnot==0){branch('D','E','F');printf("\n");}
else if(DEnot==0&&EFnot==1){branch('D','E','o');printf("\n");}
else if(DEnot==1&&EFnot==0){branch('o','E','F');printf("\n");}
col('B');printf("\n");
if(GHnot==1&&HInot==1){}
else if(GHnot==0&&HInot==0){branch('G','H','I');printf("\n");}
else if(GHnot==0&&HInot==1){branch('G','H','o');printf("\n");}
else if(GHnot==1&&HInot==0){branch('o','H','I');printf("\n");}
}

printf("**Please note that:\n any resistance value =%.1lf is considered OPEN circuit impedance**\n**any impedance input of 0 0 (format1) or 0 0 0 (format2) isconsidered short circuit impedance**\nmake anything you want!\n\n",open);
printf("I can show you some demo example about how and what you can make instructing with values(enter 1 for demo 0 to proceed.)\n");
int demomode=0;
scanf("%d",&demomode);

//****************************************************************************************************************Taking input*************************************************************
        for(int proceed=0;proceed!=1||demomode==1;){
printf("\n\nIn which format you want to input the impedances?\n\nEnter 1 for (resistance(in ohm)<space>reactance(ohm)) format or 2 for resistance(ohm)<space>inductance(millihenry)<space>capacitance(microfarad)\n");
int choice;
if (demomode==1){
    printf("\n\nDemo mode turned on.Putting values Autometically!\n");
    choice=1;
}else{scanf("%d",&choice);}
printf("***Note that****\n\nSource polarity is considered here:\n");
printf("A->B :+ve(voltage rise)\n");
printf("B->C :+ve(voltage rise)\n");
printf("C->F :+ve(voltage rise)\n");
printf("F->E :+ve(voltage rise)\n");
printf("B->E :+ve(voltage rise)\n");
printf("E->D :+ve(voltage rise)\n");
printf("D->A :+ve(voltage rise)\n");
printf("F->I :+ve(voltage rise)\n");
printf("I->H :+ve(voltage rise)\n");
printf("H->E :+ve(voltage rise)\n");
printf("H->G :+ve(voltage rise)\n");
printf("G->D :+ve(voltage rise)\n");
printf("\n*(These will be represented with sign across sources in the main diagram in the next version)*\n\n\n");
void capindform(){
    double f,mod,arg;
    printf("Input Frequency(Hz)\n");
    scanf("%lf",&f);
    for (int i=1;i<=12;i++){
    double R=0.0,L=0.0,C=0.0,imag_part=0.0;
    printf("Enter Z%d in the format of (resistance inductance capacitace):",i);
    scanf("%lf %lf %lf",&R,&L,&C);
           C *= 1e-6;
           L *= 1e-3;
           if(C==0){imag_part = 2 * M_PI * f * L ;}
           else if(L==0){imag_part =- 1 / (2 * M_PI * f * C);}
           else if(C==0&&L==0){imag_part=0;}
           else{imag_part = 2 * M_PI * f * L - 1 / (2 * M_PI * f * C);}

        Z[i] = R + imag_part * I;
    printf("Enter E%d in the format of (modulus argument):",i);
    scanf("%lf %lf",&mod,&arg);
    Z[i+12] =mod*cos(arg)+I*mod*sin(arg);
}

}
void impform(){
for (int i=1;i<=12;i++){
        double mod,arg;
    double real_part, imag_part;
        if(demomode==1){
            Z[i]=2+3*I;
            mod=2;arg=3;
            Z[3]=6113 +0*I;
            Z[7]=6113+0*I;
            Z[1]=0+0*I;
            Z[13]=0+0*I;
            Z[2]=0+0*I;
            Z[14]=0+0*I;
            Z[5]=0+0*I;
            Z[17]=0+0*I;
        }
        else{
    printf("Enter Z%d in the format of (real img):",i);
    scanf("%lf %lf",&real_part,&imag_part);
        Z[i] = real_part + imag_part * I;
    printf("Enter E%d in the format of (modulus argument):",i);
    scanf("%lf %lf",&mod,&arg);}
    sprintf(buffer, "%.2lf<%.2lf",mod,arg);
    strcpy(storedZ[i+12],buffer);
    arg=arg*(M_PI/180);
    Z[i+12] =mod* cos(arg)+I*mod*sin(arg);
}}
if(choice==1){impform();}
else{capindform();}
if(creal(Z[1])==open){ Z[2] =open+cimag( Z[2] )*I;}
if(creal(Z[2])==open){ Z[1] =open+cimag( Z[1] )*I;}
if(creal(Z[5])==open){ Z[6] =open+cimag( Z[6] )*I;}
if(creal(Z[6])==open){ Z[5] =open+cimag( Z[5] )*I;}
if(creal(Z[9])==open){ Z[10]=open+cimag( Z[10])*I;}
if(creal(Z[10])==open){Z[9] =open+cimag( Z[9] )*I;}
if(creal(Z[11])==open){Z[12]=open+cimag( Z[12])*I;}
if(creal(Z[12])==open){Z[11]=open+cimag( Z[11])*I;}

for(int i=0;i<=12;i++){if(creal(Z[i])==open){Z[i+12]=open+cimag(Z[i+12])*I;}
                       if(creal(Z[i+12])==open){Z[i-12]=open+cimag(Z[i-12])*I;}}

if(creal(Z[2])==open){ ABnot=1;vanish(1);}    else{ABnot=0;}
if(creal(Z[5])==open){ BCnot=1;vanish(2);}    else{BCnot=0;}
if(creal(Z[4])==open){ DEnot=1;combine(1,4);} else{DEnot=0;}
if(creal(Z[7])==open){ EFnot=1;combine(2,3);} else{EFnot=0;}
if(creal(Z[11])==open){GHnot=1;vanish(4);}    else{GHnot=0;}
if(creal(Z[9])==open){ HInot=1;vanish(3);}    else{HInot=0;}
if(creal(Z[3])==open){ BEnot=1;combine(1,2);} else{BEnot=0;}
if(creal(Z[8])==open){ EHnot=1;combine(3,4);} else{EHnot=0;}
// ****************************************************************************************************************************************
//finding number of equations
if(ABnot==1&&DEnot==1){GHnot=1;}if(ABnot==1&&GHnot==1){DEnot=1;}if(GHnot==1&&DEnot==1){ABnot=1;}
if(BCnot==1&&EFnot==1){HInot=1;}if(BCnot==1&&HInot==1){EFnot=1;}if(EFnot==1&&HInot==1){BCnot=1;}
if(ABnot==1&&BEnot==1){BCnot=1;}if(ABnot==1&&BCnot==1){BEnot=1;}if(BEnot==1&&BCnot==1){ABnot=1;}
if(GHnot==1&&EHnot==1){HInot=1;}if(GHnot==1&&HInot==1){EHnot=1;}if(EHnot==1&&HInot==1){GHnot=1;}

if((ABnot+BCnot+DEnot+EFnot+GHnot+HInot+BEnot+EHnot)==0) {
    unknown=4;}
else if((ABnot+BCnot+DEnot+EFnot+GHnot+HInot+BEnot+EHnot)==1){
    unknown=3;}
else if((ABnot+BCnot+DEnot+EFnot+GHnot+HInot+BEnot+EHnot)==3||(ABnot+BCnot+DEnot+EFnot+GHnot+HInot+BEnot+EHnot)==2)
   {unknown=2;}
else {unknown=1;}

for (int i=1;i<=12;i++){
    if (cimag(Z[i]) < 0) {
        sprintf(buffer, "%.2lf-j%.2lf", creal(Z[i]), -cimag(Z[i]));
            }
     else {
        sprintf(buffer, "%.2lf+j%.2lf", creal(Z[i]), cimag(Z[i]));}
    strcpy(storedZ[i],buffer);
        }

result=0;
if(demomode==1){diagram();demomode=0;
printf("Here Z1 E1 and other elements which are shorted are given input like this:0 0\nZ3 and Z7 are given:6113 0\nother input were given:2 3\n\nDemo ended.Now make using your creativity!\n");
}
else{
printf("\n\n\n\n\n\n\n\n\nIs this the circuit you want me to solve?\n\n\n\n");
diagram();
printf("Enter 1 to proceed \n");
scanf("%d",&proceed);
}
        }

printf("\n\n\n\n\n\nthis is z7 %lf\n\n\n\n\n\n",creal(Z[7]));

void putValue(){
matrix[0][0]=Z[1]+ Z[2]+ Z[3] + Z[4];  matrix[0][1]=-Z[3];                   matrix[0][2]=0+0*I;                     matrix[0][3]=-Z[4];
matrix[1][0]=-Z[3];                    matrix[1][1]=Z[3]+ Z[5]+Z[6]+Z[7];    matrix[1][2]=-Z[7];                     matrix[1][3]=0+0*I;
matrix[2][0]=0+0*I;                    matrix[2][1]=-Z[7];                   matrix[2][2]=Z[7]+ Z[8]+Z[9]+ Z[10];    matrix[2][3]=-Z[8];
matrix[3][0]=-Z[4];                    matrix[3][1]=0+0*I;                   matrix[3][2]=-Z[8];                     matrix[3][3]=Z[11] + Z[12] + Z[8] + Z[4];

consmatrix4x1[0]=Z[13] + Z[14]+ Z[15]+ Z[16];
consmatrix4x1[1]=Z[17]+ Z[18]+ Z[19]-Z[15];
consmatrix4x1[2]=Z[20]+ Z[21]+ Z[22]-Z[19];
consmatrix4x1[3]=Z[23] + Z[24]- Z[20]- Z[16];
}
putValue();
if(unknown==4||unknown==3||unknown==2){
printf("\n\n\n\n\n\nSolving with Cramer's rule will be easier for this circuit\nForming Basic Equation(matrix) for cramer's...\n\n");
        pm4(matrix);
} else {printf("\n\n\n\n\nSolving with kirchhoff's law will be easier for this circuit\n\n");}
printf("Calculating...\n\n\n\n");


void calculateCurrents() {
if (unknown==4){
    double complex det, det1, det2, det3, det4;
    double complex numMatrix[4][4];
    printf("Forming D for cramer's\n\n");
    pm4(matrix);
    printf("\n\nCalculating it's determinant by breaking down step by step...\n\n");
    det = det4x4(matrix);
    space(70);printf("Determinant of D=(%.2lf + %.2lfi)",creal(det), cimag(det));space(70);printf("\n\n");

    for (int i = 0; i < 4; i++) {
        memmove(numMatrix, matrix, sizeof(matrix));

        numMatrix[0][i] = consmatrix4x1[0];
        numMatrix[1][i] = consmatrix4x1[1];
        numMatrix[2][i] = consmatrix4x1[2];
        numMatrix[3][i] = consmatrix4x1[3];
        switch (i) {
            case 0:
                printf("Forming D1 replacing column 1\n\n");
                pm4(numMatrix);
                det1 = det4x4(numMatrix);
                space(70);printf("  Determinant of D1=(%.2lf + %.2lfi)  ",creal(det1), cimag(det1));space(70);printf("\n\n");
                break;
            case 1:
                printf("Forming D2 replacing column 2\n\n");
                pm4(numMatrix);
                det2 = det4x4(numMatrix);
                space(70);printf("  Determinant of D2=(%.2lf + %.2lfi)  ",creal(det2), cimag(det2));space(70);printf("\n\n");
                break;
            case 2:
                printf("Forming D3 replacing column 3\n\n");
                pm4(numMatrix);
                det3 = det4x4(numMatrix);
                space(70);printf("  Determinant of D3=(%.2lf + %.2lfi)  ",creal(det3), cimag(det3));space(70);printf("\n\n");
                break;

            case 3:
                printf("Forming D4 replacing column 4\n\n");
                pm4(numMatrix);
                det4 = det4x4(numMatrix);
                space(70);printf("  Determinant of D4=(%.2lf + %.2lfi)  ",creal(det4), cimag(det4));space(70);printf("\n\n");
                break;}
                            }
                    printf("Calculating Currents[Ix]=[Dx/D]...\n\n");
    curren[1]= det1/det;
    curren[5] =det2/det;
    curren[9] =det3/det;
    curren[12]=det4/det;
    }
            else if(unknown==3){

                int zno,zno2,con1,con2,con3,i1,i2,i3,i4,ic;
                double complex det, det1, det2, det3;
                double complex numMatrix[3][3];
                double complex constant[10];

                    if (GHnot==1){zno=11;zno2=12;con1=1;con2=2;con3=3;i1=1;i2=5;i3=9;i4=12;}
               else if (HInot==1){zno=9;zno2=10;con1=1;con2=2;con3=4;i1=1;i2=5;i3=12;i4=9;}
               else if (BCnot==1){zno=5;zno2=6;con1=1;con2=3;con3=4;i1=1;i2=9;i3=12;i4=5;}
               else if (ABnot==1){zno=2;zno2=1;con1=2;con2=3;con3=4;i1=5;i2=9;i3=12;i4=1;}
               else if (EHnot==1){zno=8;zno2=8;con1=1;con2=2;con3=8;i1=1;i2=5;i3=9;i4=12;ic=9;}
               else if (BEnot==1){zno=3;zno2=3;con1=5;con2=3;con3=4;i1=1;i2=9;i3=12;i4=5;ic=1;}
               else if (DEnot==1){zno=4;zno2=4;con1=6;con2=2;con3=3;i1=1;i2=5;i3=9;i4=12;ic=1;}
               else if (EFnot==1){zno=7;zno2=7;con1=1;con2=7;con3=4;i1=1;i2=5;i3=12;i4=9;ic=5;}

                              Z[zno]=0+0*I;
                              Z[zno+12]=0+0*I;
                              Z[zno2]=0+0*I;
                              Z[zno2+12]=0+0*I;
                              putValue();

                              if (vanEqn[0]!=0){
                            printf("Forming matrix reducing col and row no:%d ...\n\n",vanEqn[0]);
                                pc(consmatrix4x1[0]);printf("\n");pc(consmatrix4x1[1]);printf("\n");pc(consmatrix4x1[2]);printf("\n");pc(consmatrix4x1[3]);printf("\n");printf("\n");printf("\n");
                            make3x3(vanEqn[0]);     make3x1(vanEqn[0]);
                                                                    pc(consmatrix3x1[0]);printf("\n");pc(consmatrix3x1[1]);printf("\n");pc(consmatrix3x1[2]);printf("\n");
                    }
                    else if (comEqn1[0]!=0){
                            AddCol4x4(comEqn1[0],comEqn2[0]);    AddRow4x(comEqn1[0],comEqn2[0]);
                                                                             make3x1(comEqn2[0]);
                            printf("Forming matrix modifying column...\n\n");
                            pm4(matrix);
                            AddRow4x4(comEqn1[0],comEqn2[0]);
                            printf("Forming matrix modifying Row...\n\n");
                            pm4(matrix);
                            make3x3(comEqn2[0]);}

                printf("Forming matrix from the equation...\n\n");
                pm3(reduced);
                det = det3x3(reduced);
                printf("\nDeterminant of original matrix (D)\n");
                pc(det);

    for (int h = 0; h < 3; h++) {
        memmove(numMatrix,reduced,sizeof(reduced));
        printf("\nForming matrix by replacing column no. %d with constants for calculating D(%d)...\n\n",h+1,h+1);
        numMatrix[0][h]=consmatrix3x1[0];
        numMatrix[1][h]=consmatrix3x1[1];
        numMatrix[2][h]=consmatrix3x1[2];

        pm3(numMatrix);
        switch (h) {
            case 0:
                det1 = det3x3(numMatrix);
                space(70);printf("Calculated Determinant of this 3x3 matrix(D%d): (%.2lf + %.2lfi)\n",h+1,creal(det1), cimag(det1));
                break;
            case 1:
                det2 = det3x3(numMatrix);
                space(70);printf("Calculated Determinant of this 3x3 matrix(D%d): (%.2lf + %.2lfi)\n",h+1,creal(det2),cimag(det2));
                break;
            case 2:
                det3 = det3x3(numMatrix);
                space(70);printf("Calculated Determinant of this 3x3 matrix(D%d): (%.2lf + %.2lfi)\n",h+1,creal(det3), cimag(det3));
                break;}
    }
            curren[i1]=det1/det;printf("\n\nCurrent through Z%d=(D1)/D                     calculated!\n\n",i1);
            pc(curren[i1]);
            curren[i2]=det2/det;printf("\n\nCurrent through Z%d=(D2)/D                     calculated!\n\n",i2);
            pc(curren[i2]);
            curren[i3]=det3/det;printf("\n\nCurrent through Z%d=(D3)/D                     calculated!\n\n",i3);
            pc(curren[i3]);printf("\n\n");
            if(comEqn1[0]!=0){curren[i4]=curren[ic];}
                    else{ curren[i4]= 0+0*I;}

           Z[zno]= open;
           Z[zno+12]=open;
           Z[zno2]=open;
           Z[zno2+12]=open;
                }
        else if (unknown==2){
                double complex numerator[2][2];
                double complex det,det1,det2;
                int combination=0,current1,current2;
                int indicator[30];
                 for(int i=1;i<=24;i++){
                            indicator[i]=0;
                            if(creal(Z[i])==open){Z[i]=0+0*I;
                                                indicator[i]=i;}
                        }
                        putValue();

        if (vanEqn[0]!=0){
                if(vanEqn[1]!=0){

                    make3x3(vanEqn[0]);        make3x1(vanEqn[0]);
                    if (vanEqn[1]>vanEqn[0]){vanEqn[1]--;}
                    make2x2(vanEqn[1]);     make2x1(vanEqn[1]);}

                else if (vanEqn[0]==comEqn1[0]){
                                make3x3(vanEqn[0]);      make3x1(vanEqn[0]);
                                if (comEqn1[0]>vanEqn[0]){comEqn1[0]--;}
                                make2x2(comEqn1[0]);  make2x1(comEqn1[0]);
                            }

                else if(vanEqn[0]==comEqn2[0]){
                                make3x3(vanEqn[0]);          make3x1(vanEqn[0]);
                                if (comEqn2[0]>vanEqn[0]){comEqn2[0]--;}
                                make2x2(comEqn2[0]);      make2x1(comEqn2[0]);}
                else {  combination=1;
                            AddCol4x4(comEqn1[0],comEqn2[0]);
                            printf("Forming matrix modifying column...\n\n");
                            pm4(matrix);
                            AddRow4x4(comEqn1[0],comEqn2[0]);       AddRow4x(comEqn1[0],comEqn2[0]);
                            printf("Forming matrix modifying Row...\n\n");
                            pm4(matrix);
                            make3x3(comEqn2[0]);                           make3x1(comEqn2[0]);
                            if (vanEqn[0]>comEqn2[0]){vanEqn[0]--;}
                            make2x2(vanEqn[0]);                         make2x1(vanEqn[0]);
                        }}

                else{    combination=2;
                            AddCol4x4(comEqn1[0],comEqn2[0]);
                            printf("Forming matrix modifying column...\n\n");
                            pm4(matrix);
                            AddRow4x4(comEqn1[0],comEqn2[0]);        AddRow4x(comEqn1[0],comEqn2[0]);
                            printf("Forming matrix modifying Row...\n\n");
                            pm4(matrix);

                            AddCol4x4(comEqn1[1],comEqn2[1]);   AddRow4x(comEqn1[1],comEqn2[1]);
                            printf("Forming matrix modifying column...\n\n");
                            pm4(matrix);
                            AddRow4x4(comEqn1[1],comEqn2[1]);
                            printf("Forming matrix modifying Row...\n\n");
                            pm4(matrix);

                            make3x3(comEqn2[0]);   make3x1(comEqn2[0]);
                            if (comEqn2[1]>comEqn2[0]){comEqn2[1]--;}
                            make2x2(comEqn2[1]);          make2x1(comEqn2[1]);
                }
                det=det2x2(decreased[0][0],decreased[0][1],decreased[1][0],decreased[1][1]);

             for (int h = 0; h < 2; h++) {
        memmove(numerator,decreased,sizeof(decreased));
        printf("\n\n\nForming matrix by replacing column no. %d with constants for calculating D(%d)...\n\n",h+1,h+1);
        numerator[0][h]=consmatrix2x1[0];
        numerator[1][h]=consmatrix2x1[1];
        pm2(numerator);

                switch (h) {
            case 0:
                det1 = det2x2(numerator[0][0],numerator[0][1],numerator[1][0],numerator[1][1]);
                printf("\n\nCalculated Determinant of this 2x2 matrix(D%d): (%.2lf + %.2lfi)\n",h+1,creal(det1), cimag(det1));
                break;
            case 1:
                det2 = det2x2(numerator[0][0],numerator[0][1],numerator[1][0],numerator[1][1]);
                printf("\n\nCalculated Determinant of this 3x3 matrix(D%d): (%.2lf + %.2lfi)\n",h+1,creal(det2), cimag(det2));
                break;}
    }
                        int a=0,b=0,c=0,d=0,x=0,y=0;

                        if(combination==0){
                            a=vanEqn[0];
                            b=vanEqn[1];
                            for (int i=0;i<=4;i++){
                                if(i!=a&&i!=b){
                                    if (c==0){c=i;}
                                    else if (d==0){d=i;}
                                }
                            }
                            arrange(c,d,0);
                            iloop[highest]=det1/det;
                            iloop[medium]=det2/det;
                        }
                        else if (combination==1){
                            a=vanEqn[0];
                            b=comEqn1[0];
                            c=comEqn2[0];
                            for(int i=0;i<=4;i++){
                                if ((i!=vanEqn[0])&&(i!=comEqn1[0])&&(i!=comEqn2[0])){
                                    d=i;
                                }
                                arrange(b,d,0);
                                iloop[highest]=det1/det;
                                iloop[medium]=det2/det;
                                iloop[c]=iloop[b];
                            }
                        }
                        else if(combination==2){
                                a=comEqn1[0];b=comEqn2[0];c=comEqn1[1];d=comEqn2[1];
                            if(a==c||a==d||b==c||b==d){
                                x=a;
                                for(int i=1;i<=4;i++){
                                    if(i!=a&&i!=b&&i!=c&&i!=d){
                                        y=i;
                                    }
                                }
                                current1=x;current2=y;
                                arrange(current1,current2,0);
                                iloop[highest]=det1/det;
                                iloop[medium]=det2/det;
                                for (int j=1;j<=4;j++){
                                    if(j!=y){iloop[j]=iloop[current1];}
                                }}

                        else{current1=a;current2=c;}
                            arrange(current1,current2,0);
                            iloop[highest]=det1/det;
                            iloop[medium]=det2/det;
                            iloop[b]=iloop[a];
                            iloop[d]=iloop[c];
                        }
                        curren[1]= iloop[1];
                        curren[5]= iloop[2];
                        curren[9]= iloop[3];
                       curren[12]= iloop[4];
                       for(int i=0;i<=24;i++){
                        if (indicator[i]!=0){
                            Z[i]=open;}}
                    printf("\n\nCalculating currents[D(x)/D]...\n\n");
                            }

        else {
            double complex totalZ,totalE;
            totalZ=0+0*I;
            totalE=0+0*I;
            int indicator[30];
            for(int t=1;t<=12;t++){ indicator[t]=0;indicator[t+12]=0;
                    if(creal(Z[t])==open){Z[t]=0+ 0*I;indicator[t]=t;}
                    if(creal(Z[t+12])==open){Z[t+12]=0+0*I;indicator[t+12]=t+12;}
                    putValue();
                totalZ +=Z[t];
                totalE +=Z[t+12];}
                curren[1]=totalE/totalZ;
                for (int i=1;i<=12;i++){
                    curren[i]=curren[1];
                    if(indicator[i]!=0){Z[i]=open;}
                    if(indicator[i+12]!=0){Z[i+12]=open;}
                }

                printf("\nTotal voltage rise =(%.2lf + %.2lf)\n",creal(totalE),cimag(totalE));
                printf("\nEquivalant impedance =(%.2lf + %.2lf)\n\n",creal(totalZ),cimag(totalZ));
                printf("Calculating currents[T.voltage rise/Equivalant Impedance]...");

                }
}
calculateCurrents();
void putcurrents(){
curren[2]= curren[1];            curren[6]=  curren[5];              curren[10]=curren[9];             curren[11]=curren[12];
curren[3]= curren[5]-curren[1];  curren[8]=  curren[9]-curren[12];   curren[7]=curren[9]-curren[5];    curren[4]=curren[12]-curren[1];
curren[9]=-curren[9];            curren[10]=-curren[10];             curren[6]=-curren[6];             curren[11]=-curren[11];}

if(unknown!=1){putcurrents();}
for(int r=1;r<=12;r++){
    V[r]=Z[r]*curren[r];
                }
double complex poler(double complex rect){
    double r = cabs(rect);
    double theta =(180.0 / M_PI)*carg(rect);
    return r+theta*I;
}
for(int i=1;i<=12;i++){
        curren[i]=poler(curren[i]);
        V[i]=poler(V[i]);
}

for (int i = 1; i <= 12; i++) {
        printf("Current I[%d] = %.2lf <%.2lf\n", i, creal(curren[i]), cimag(curren[i]));
        printf("Voltage V[%d]: %.2lf < %.2lf\n", i, creal(V[i]), cimag(V[i]));
                }
for (int i=1;i<=12;i++){
   sprintf(storedV[i],"%.2lf<%.2lf",creal(V[i]),cimag(V[i]));
   sprintf(storedV[i+12],"%.2lf<%.2lf",creal(curren[i]),cimag(curren[i]));}

result=1;
printf("\n\n\n\n\n\nHere is your circuit with results!\n\n\n\n");
diagram();
printf("***Note that****\n\nSource polarity is considered here:\n");
printf("A->B :+ve(voltage rise)\n");
printf("B->C :+ve(voltage rise)\n");
printf("C->F :+ve(voltage rise)\n");
printf("F->E :+ve(voltage rise)\n");
printf("B->E :+ve(voltage rise)\n");
printf("E->D :+ve(voltage rise)\n");
printf("D->A :+ve(voltage rise)\n");
printf("F->I :+ve(voltage rise)\n");
printf("I->H :+ve(voltage rise)\n");
printf("H->E :+ve(voltage rise)\n");
printf("H->G :+ve(voltage rise)\n");
printf("G->D :+ve(voltage rise)\n");
printf("\n*(These will be represented with sign across sources in the main diagram in the next version)*\n\n");

printf("\n\n\nThank You!\nDeveloper : Ashraful Haque Jisan\n            EEE KUET,BATCH 2K22\n\n");
printf("\nThis software and associated files are the intellectual property of the copyright holder.\n");
printf("Unauthorized copying, modification, distribution, or use of this software or its components, in any form, is strictly prohibited without prior written permission from the copyright holder.\n\n\n");
printf("Please mail at ashrafulhaquejisan@gmail.com for any feedback.\n\n");
double exit;
printf("ENTER YOUR FAVOURITE NUMBER TO EXIT   \n\n");
for(int f=0;exit!=112.0;f++){
    scanf("%lf",&exit);
    if(exit!=112.0){
        printf("Wrong Answer        HINT!(Your favourite number should be 112)  ");
    }
}
  return 0;
}
         void space(int s) {
            for (int i=0;i<s;i++){
            printf(" ");}}

         void curl(int c){
            for(int i=0;i<c;i++){
            printf("~");}}

         void dash(int d){
            for(int i=0;i<d;i++){
            printf("-");}}

         void star(int t,char a,char b){
            printf("%c",a);
            for(int j=0;j<t;j++){
                    printf("*");}
            printf("%c",b);}

        void volt(int n){
                int len =strlen(storedV[n]);
                int totalSpaces=13-len;
                int halfSpaces =totalSpaces/2;

                for (int i=0;i<halfSpaces;i++) {
                            printf(" ");}
                            printf("%s", storedV[n]);
                        for (int i = 0; i < (totalSpaces - halfSpaces); i++) {
                            printf(" ");}}

        void curr(int n){volt(n+12);}

        void imp(int t,int n){
                int len = strlen(storedZ[n]);
                int totalSpaces =14-len;
                int halfSpaces = totalSpaces / 2;

                    if(t==1){
                            if(creal(Z[n])==0&&cimag(Z[n])==0){star(14,'*','*');}
                            else{printf("|");

                        for (int i=0;i<halfSpaces;i++) {
                            printf(" ");}
                            printf("%s", storedZ[n]);
                        for (int i = 0; i < (totalSpaces - halfSpaces); i++) {
                            printf(" ");}
                            printf("|");}
                            cellcall=0;}

                    else if(t==2){int o,p;

                        void reptop(int type){
                                        if (n==1){o=3;p=6;}else if(n==12){o=8;p=10;}else if (n==13){o=15;p=18;}else if(n==24){o=20;p=22;} else{}
                            void part(int j){ int w=53,e=69,t=50;if(j==p){w=0;e=0;t=0;}
                                        if(cabs(Z[j])==0){space(15);printf("*");space(w);} else if(creal(Z[j])==open){space(e);}
                                            else {if(result==1&&cellcall!=1){space(13);if(type==1){printf("--");}else{printf("| ");}}else {space(15);}if(type==1){curl(4);}else{printf("|  |");}space(t);}
                                        }
                                part(n);
                                part(o);
                                part(p);printf("\n");}

                        void repmid(int z){int c=53,a=69,t=39;if(z==p){c=0;a=0;t=0;}
                        if(cabs(Z[z])==0){space(15);printf("*");space(c);} else if(creal(Z[z])==open){space(a);} else{if(result==1&&cellcall!=1){volt(z);printf("V ");}else {space(15);}
                                                                                            len = strlen(storedZ[z]);
                                                                                            totalSpaces =14-len;
                                                                                            printf("|");
                                                                                            printf("%s", storedZ[z]);
                                                                                            for (int i = 0; i < totalSpaces; i++) {
                                                                                            printf(" ");}space(t);
                                                                                            }}
                                        reptop(1);reptop(2);repmid(n);repmid(o);repmid(p);printf("\n");
                                        reptop(2);reptop(1);
                            cellcall=0;}}

           void cell(int t,int n) {cellcall=1;imp(t,n+12);}

           void branch(char a,char b,char c){
                void roof(int f){printf("*");space(13);if(isImp1!=1||isImp2!=1){space(14);}else{curl(14);}space(14);if(isCell1!=1||isCell2!=1){space(14);}else{dash(14);}space(13);}
                 void bottom(int b){printf("*");space(11);if(isImp1!=1||isImp2!=1){space(30);}else{if(result==1){printf("| ");}
                                        else {space(2);}curl(14);if(result==1){printf(" |");}else {space(2);}space(12);}
                        if(isCell1!=1||isCell2!=1){space(27);}else{dash(14);space(13);}}
                void rating(int r){
                    printf("*");space(11);if(result==1){if(isImp1!=1||isImp2!=1){space(18);}
                                            else{printf("|_");volt(r);printf("V");printf("_|");}space(3);dash(14);printf(">");curr(r);printf("amp");space(5);}else{space(57);}}
                int z,j;
                char t,s;

                if(c=='o'||a=='o'){
                    if(c=='o'){
                            j=15;t=a;s=b;
                            if(a=='A'){z=2;}
                            else if (a=='D'){z=4;}
                            else if (a=='G'){z=11;}
                            }
                    if(a=='o'){
                        j=84;t=b;s=c;
                    if(b=='B'){z=5;}
                    else if(b=='E'){z=7;}
                    else if(b=='H'){z=9;}
                            }
                                    if(cabs(Z[z])==0){isImp1=0;}else{isImp1=1;isImp2=1;}if(cabs(Z[z+12])==0){isCell1=0;}else{isCell1=1;isCell2=1;}
                                       if(j==84){space(15);printf("*");j=j-16;}space(j);roof(1);printf("*");if(j==15){space(68);printf("*\n");}else{printf("\n");}
                                       if(j==68){space(15);printf("*");}space(j);star(1,t,'*');star(8,'*','*');imp(1,z);star(10,'*','*');cell(1,z);star(10,'*','*');printf("%c",s);
                                                                                            if(j==15){space(68);printf("*\n");}else{printf("\n");}
                                       if(j==68){space(15);printf("*");}space(j);bottom(11);printf("*");if(j==15){space(68);printf("*\n");}else{printf("\n");}
                                       if(j==68){space(15);printf("*");}space(j);rating(z);printf("*");if(j==15){space(68);printf("*");}else{}}

                else{int z,x,i1,i2;

                            if(a=='A'){z=2;x=5;}
                            else if (a=='D'){z=4;x=7;}
                            else if (a=='G'){z=11;x=9;}
                            if(cabs(Z[z])==0){imp1=0;}else{imp1=1;}if(cabs(Z[z+12])==0){cell1=0;}else{cell1=1;}
                            if(cabs(Z[x])==0){imp2=0;}else{imp2=1;}if(cabs(Z[x+12])==0){cell2=0;}else{cell2=1;}

                            space(15);isImp1=imp1;isCell1=cell1;roof(1);isImp1=1;isCell1=1;
                            isImp2=imp2;isCell2=cell2;roof(1);isImp2=1;isCell2=1;printf("*\n");space(15);

                                    star(1,a,'*');star(8,'*','*');imp(1,z);star(10,'*','*');cell(1,z);star(10,'*','*');
                                    star(1,b,'*');star(8,'*','*');imp(1,x);star(10,'*','*');cell(1,x);star(10,'*','*');printf("%c\n",c);
                                    space(15);    isImp1=imp1;isCell1=cell1;bottom(1);isImp1=1;isCell1=1;
                                                  isImp2=imp2;isCell2=cell2;bottom(1);isImp2=1;isCell2=1;printf("*\n");
                                    space(15);isImp1=imp1;isCell1=cell1;rating(z);isImp1=1;isCell1=1;
                                              isImp2=imp2;isCell2=cell2;rating(x);isImp2=1;isCell2=1;printf("*");}}
            void col(char a){
                            int c,d,e;
                                if (a=='A'){c=1;d=3;e=6;
                                imp(2,1);}else{imp(2,12);c=12;d=8;e=10;}

                                                space(15);if (creal(Z[c])==open){space(69);}else{
                                                printf("*");if(result==1){space(5);printf("^");space(62);}else{space(68);}}
                                         if (creal(Z[d])==open){space(69);}else{
                                                printf("*");if(result==1){space(5);printf("^");space(62);}else{space(68);}}
                                         if (creal(Z[e])==open){printf("\n");}else{if(result==1){printf("*    ^\n");}else{printf("*\n");}}

                                space(15);if (creal(Z[c])==open){space(69);}else{
                                                printf("*");if(result==1){space(5);printf("|");curr(c);printf("amp");space(46);}else{space(68);}}
                                          if (creal(Z[d])==open){space(69);}else{
                                          printf("*");if(result==1){space(5);printf("|");curr(d);printf("amp");space(46);}else{space(68);}}
                                          if (creal(Z[e])==open){printf("\n");}else{if(result==1){printf("*    |");curr(e);printf("amp\n");}else{printf("*\n");}}

                                if (a=='A'){
                                cell(2,1);}else{cell(2,12);}
                                space(15);if (creal(Z[c])==open){space(69);}else{printf("*");space(68);}if (creal(Z[d])==open){space(69);}else{printf("*");space(68);}
                                if (creal(Z[e])==open){}else{printf("*");}}
