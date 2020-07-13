float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}

void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

void calculateGammaA(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float x1, x2, x3, x4;
    x1 = n1.getY();
    x2 = n2.getY();
    x3 = n3.getY();
    x4 = n4.getY();

	G1.at(0).at(0) = -(3*x1-x2-x3-x4+3*z1-z2-z3-z4-5);   	G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = -(4*x1-2*x2-x3-x4+4*z1-2*z2-z3-z4-5);  G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = -(4*x1-x2-2*x3-x4+4*z1-z2-2*z3-z4-5);  G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = -(4*x1-x2-x3-2*x4+4*z1-z2-z3-2*z4-5);  G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = -(3*x1-x2-x3-x4+3*z1-z2-z3-z4-5);   	G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = -(4*x1-2*x2-x3-x4+4*z1-2*z2-z3-z4-5);  G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = -(4*x1-x2-2*x3-x4+4*z1-z2-2*z3-z4-5); 	G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = -(4*x1-x2-x3-2*x4+4*z1-z2-z3-2*z4-5);  G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = -(3*x1-x2-x3-x4+3*z1-z2-z3-z4-5); 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = -(4*x1-2*x2-x3-x4+4*z1-2*z2-z3-z4-5); 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = -(4*x1-x2-2*x3-x4+4*z1-z2-2*z3-z4-5);
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = -(4*x1-x2-x3-2*x4+4*z1-z2-z3-2*z4-5); 
}

void calculateGammaC(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float x1, x2, x3, x4;
    x1 = n1.getY();
    x2 = n2.getY();
    x3 = n3.getY();
    x4 = n4.getY();

	G1.at(0).at(0) = -(6*x1-2*x2-2*x3-2*x4+3*z1-z2-z3-z4-25);  	G1.at(0).at(1) = 0;   										G1.at(0).at(2) = 0;
	G1.at(1).at(0) = -(8*x1-4*x2-2*x3-2*x4+4*z1-2*z2-z3-z4-25); G1.at(1).at(1) = 0;   										G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = -(8*x1-2*x2-4*x3-2*x4+4*z1-z2-2*z3-z4-25); G1.at(2).at(1) = 0;   										G1.at(2).at(2) = 0;
	G1.at(3).at(0) = -(8*x1-2*x2-2*x3-4*x4+4*z1-z2-z3-2*z4-25); G1.at(3).at(1) = 0;   										G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   										G1.at(4).at(1) = -(6*x1-2*x2-2*x3-2*x4+3*z1-z2-z3-z4-25);   G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   										G1.at(5).at(1) = -(8*x1-4*x2-2*x3-2*x4+4*z1-2*z2-z3-z4-25); G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   										G1.at(6).at(1) = -(8*x1-2*x2-4*x3-2*x4+4*z1-z2-2*z3-z4-25); G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   										G1.at(7).at(1) = -(8*x1-2*x2-2*x3-4*x4+4*z1-z2-z3-2*z4-25); G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   										G1.at(8).at(1) = 0;   										G1.at(8).at(2) = -(6*x1-2*x2-2*x3-2*x4+3*z1-z2-z3-z4-25);   
	G1.at(9).at(0) = 0;   										G1.at(9).at(1) = 0;   										G1.at(9).at(2) = -(8*x1-4*x2-2*x3-2*x4+4*z1-2*z2-z3-z4-25);
    G1.at(10).at(0) = 0; 										G1.at(10).at(1) = 0; 										G1.at(10).at(2) = -(8*x1-2*x2-4*x3-2*x4+4*z1-z2-2*z3-z4-25);
	G1.at(11).at(0) = 0;  										G1.at(11).at(1) = 0;  										G1.at(11).at(2) = -(8*x1-2*x2-2*x3-4*x4+4*z1-z2-z3-2*z4-25);  
}

void calculateGammaG(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float x1, x2, x3, x4;
    x1 = n1.getY();
    x2 = n2.getY();
    x3 = n3.getY();
    x4 = n4.getY();

	G1.at(0).at(0) = -(6*x1-2*x2-2*x3-2*x4+3*z1-z2-z3-z4-25);  	G1.at(0).at(1) = 0;   										G1.at(0).at(2) = 0;
	G1.at(1).at(0) = -(8*x1-4*x2-2*x3-2*x4+4*z1-2*z2-z3-z4-25); G1.at(1).at(1) = 0;   										G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = -(8*x1-2*x2-4*x3-2*x4+4*z1-z2-2*z3-z4-25); G1.at(2).at(1) = 0;   										G1.at(2).at(2) = 0;
	G1.at(3).at(0) = -(8*x1-2*x2-2*x3-4*x4+4*z1-z2-z3-2*z4-25); G1.at(3).at(1) = 0;   										G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   										G1.at(4).at(1) = -(6*x1-2*x2-2*x3-2*x4+3*z1-z2-z3-z4-25);   G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   										G1.at(5).at(1) = -(8*x1-4*x2-2*x3-2*x4+4*z1-2*z2-z3-z4-25); G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   										G1.at(6).at(1) = -(8*x1-2*x2-4*x3-2*x4+4*z1-z2-2*z3-z4-25); G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   										G1.at(7).at(1) = -(8*x1-2*x2-2*x3-4*x4+4*z1-z2-z3-2*z4-25); G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   										G1.at(8).at(1) = 0;   										G1.at(8).at(2) = -(6*x1-2*x2-2*x3-2*x4+3*z1-z2-z3-z4-25);   
	G1.at(9).at(0) = 0;   										G1.at(9).at(1) = 0;   										G1.at(9).at(2) = -(8*x1-4*x2-2*x3-2*x4+4*z1-2*z2-z3-z4-25);
    G1.at(10).at(0) = 0; 										G1.at(10).at(1) = 0; 										G1.at(10).at(2) = -(8*x1-2*x2-4*x3-2*x4+4*z1-z2-2*z3-z4-25);
	G1.at(11).at(0) = 0;  										G1.at(11).at(1) = 0;  										G1.at(11).at(2) = -(8*x1-2*x2-2*x3-4*x4+4*z1-z2-z3-2*z4-25);  
}

void calculateGammaE(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float x1, x2, x3, x4;
    x1 = n1.getY();
    x2 = n2.getY();
    x3 = n3.getY();
    x4 = n4.getY();

	G1.at(0).at(0) = -(3*z1-z2-z3-z4);  						G1.at(0).at(1) = 0;   										G1.at(0).at(2) = 0;
	G1.at(1).at(0) = -(4*z1-2*z2-z3-z4); 						G1.at(1).at(1) = 0;   										G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = -(4*z1-z2-2*z3-z4); 						G1.at(2).at(1) = 0;   										G1.at(2).at(2) = 0;
	G1.at(3).at(0) = -(4*z1-z2-z3-2*z4); 						G1.at(3).at(1) = 0;   										G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   										G1.at(4).at(1) = -(3*z1-z2-z3-z4);    						G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   										G1.at(5).at(1) = -(4*z1-2*z2-z3-z4); 						G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   										G1.at(6).at(1) = -(4*z1-z2-2*z3-z4);  						G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   										G1.at(7).at(1) = -(4*z1-z2-z3-2*z4);  						G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   										G1.at(8).at(1) = 0;   										G1.at(8).at(2) = -(3*z1-z2-z3-z4);  
	G1.at(9).at(0) = 0;   										G1.at(9).at(1) = 0;   										G1.at(9).at(2) = -(4*z1-2*z2-z3-z4); 
    G1.at(10).at(0) = 0; 										G1.at(10).at(1) = 0; 										G1.at(10).at(2) = -(4*z1-z2-2*z3-z4); 
	G1.at(11).at(0) = 0;  										G1.at(11).at(1) = 0;  										G1.at(11).at(2) = -(4*z1-z2-z3-2*z4); 
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixK,matrixG,matrixD,matrixC,matrixE,matrixW;
    float u_bar,nu,rho,Ve,J,Determinant;
    
    /* [ A+K  G ]
       [  D   0 ]
    */
    
    /* [ A+K  C ]
       [  G   E-W ]
    */

    //Matrix A
    Matrix g_matrix, Alpha, Beta, gA_matrix,gC_matrix,gG_matrix;

    u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(120*Determinant);
    //cout<<"EQUIS DEEE";
    
    calculateGammaA(gA_matrix,m,e);
    calculateGammaC(gC_matrix,m,e);
    calculateGammaG(gG_matrix,m,e);
    calculateGamma(g_matrix);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_a, productMatrixMatrix(gA_matrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);


    //Matrix K
    Matrix Alpha_t,Beta_t;

    nu = m.getParameter(DYNAMIC_VISCOSITY);
    Ve = calculateLocalVolume(e,m);
    
    float real_k = (float) (1)/(6*Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);

	//Matrix C
	Matrix Omega;
	calculateOmega(Omega);
	float real_c = (float) (J)/(24*Determinant);
	 
	productRealMatrix(real_c, productMatrixMatrix(gC_matrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixC);
    
    //Matrix G
    
    rho = m.getParameter(DENSITY);
    float real_g = (float) (J/(24*Determinant));

    
    productRealMatrix(real_g, productMatrixMatrix(gG_matrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixG);

    //Matrix D
    Matrix g_matrix_t,Omega_t;
    float real_d = (float)(J/(24*Determinant));

    transpose(Omega, Omega_t);
    transpose(g_matrix,g_matrix_t);
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);
    
    //matrix E
    
    Matrix gG_matrix_t;
    float real_e = (float)(J/(24*Determinant));

    //transpose(Omega, Omega_t);
    transpose(gG_matrix,gG_matrix_t);
    productRealMatrix(real_e,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,gG_matrix_t,3,3,12),12,3,12),matrixE);
    
    //matrix w
    float real_w = (float) (1)/(6*Determinant*Determinant);
    productRealMatrix(real_w,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixW);

    //Matrix M
    Matrix M;
    zeroes(M,24);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,23,matrixC);
    ubicarSubMatriz(M,12,23,0,11,matrixG);
    ubicarSubMatriz(M,12,23,12,23,restarMatrix(matrixE,matrixW,12,12));

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix;

    calculateF(f, m);

    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,16);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    
    return b;
}
