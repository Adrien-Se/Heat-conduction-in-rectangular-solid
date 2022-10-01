/*
Eliot CHRISTON
28605232
L3 DM Méca/EEA

TP2 - mini-projet : Conduction de la chaleur dans une plaque rectangulaire
    Ce TP fait l’objet d’un compte-rendu. Le compte-rendu a remettre sur Moodle pour le 22/04/2022, limite
    a 6 pages maximum (hors listing du programme et annexes eventuelles) doit comprendre une introduction, une
    reponse redigee aux diverses questions, le trace des figures demandees, une analyse des courbes et des resultats
    obtenus, et une conclusion generale. Les fichiers du ou des programmes devront egalement être fournis.
*/

//  Pour mieux s'y retrouver :
//      /*AFF*/ : ligne d'affichage sur la console
//      /*MAT*/ : ligne qui concerne le fichier_matlab
//      /*FI0*/ : ligne qui concerne le fichier0


//////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 40000
#define NMAX 300
#define AFF 200         // nx*ny maximum pour affichage sur console
#define NB_MAILLAGES 3  // Nombre de maillages differents   |ATTENTION| doit etre inferieur a la taille des tableaux nxs et nys


//////////////////////////////////////////////////////////////////////////////////


float norme2_vect(float x[MAX], int n) {
    float var = 0;
    int i;
    for (i = 0; i < n; i++) var += (x[i] * x[i]);
    return sqrt(var);
}

void ecriture (float Q[MAX], float T[NMAX][NMAX], float x[NMAX], float y[NMAX], int nx, int ny, float h, float gs, float gw, float ge, float gn ) {

    // VARIABLES LOCALES
    int i, j, IJ, I=0, J=0, nn = nx*ny;

    // CHANGEMENT D'INDICES
    for (IJ = 0; IJ<nn ; IJ++){
        I = (IJ%nx)+1;
        if (IJ%nx==0){
            J += 1;
        }
        T[I][J] = Q[IJ+nx]; // on resonne en termes de lignes et colonnes, I lignes, J colonnes
    }

    // CONDITIONS AUX LIMITES
    for (i=0; i<nx+2 ; i++){T[i][0]     = gs;}
    for (i=0; i<nx+2 ; i++){T[i][ny+1]  = gn;}
    for (j=1; j<ny+1 ; j++){T[0][j]     = gw;}
    for (j=1; j<ny+1 ; j++){T[nx+1][j]  = ge;}

    // AFFICHAGE CONSOLE
    /*AFF*/if (nx*ny < AFF) { //si maillage avec peu de points
    /*AFF*/    printf("\n\nT(xi,yi)");
    /*AFF*/    printf("\n----");for (i = 0; (i<nx+2) ;i++){printf("--------");}printf("\n");
    /*AFF*/    for (j = ny + 1; j >= 0 ; --j) {
    /*AFF*/        for (i = 0; i <= nx + 1; i++) {
    /*AFF*/            if (i==0){printf("%2.2f|\t", y[j]);}
    /*AFF*/            printf("%2.2f\t", T[i][j]);
    /*AFF*/        }
    /*AFF*/        printf("\n");
    /*AFF*/    }
    /*AFF*/    printf("----");  for (i = 0; (i<nx+2) ;i++){printf("--------");}printf("\n\t");
    /*AFF*/    for (i=0; i<nx+2 ;i++) printf("%2.2f\t", x[i]);
    /*AFF*/    printf("\n");
    /*AFF*/}

}


//////////////////////////////////////////////////////////////////////////////////



int main() {

    printf("\nDEBUT DU PROGRAMME\n");

/// CREATION FICHIERS OUTPUT
    /*FI0*/FILE *fichier0 = fopen("res_sor.dat","wt");  // fichier impose par l'enonce. Peut permettre de comparer des valeurs ou de les tracer avec gnuplot eventuellement. (Matlab recommande)
    /*FI0*/if(fichier0==NULL) {printf("\nPb ouverture fichier0.\n");exit(0);}

    /*MAT*/FILE *fichier_matlab = fopen("data_TP2.m","wt");
    /*MAT*/if(fichier_matlab==NULL) {printf("\nPb ouverture fichier_matlab.\n");exit(0);}

    // Initialisation du fichier_matlab qui est une fonction sous matlab
    /*MAT*/fprintf(fichier_matlab,"function [Norm,x,y,T] = data_TP2()\n\n");
    /*MAT*/fprintf(fichier_matlab,"\n%%%% INITIALISATION\n\tNorm = cell(%d,1);\n\tx = cell(%d,1);\n\ty = cell(%d,1);\n\tT = cell(%d,1);\n\n",NB_MAILLAGES,NB_MAILLAGES,NB_MAILLAGES,NB_MAILLAGES);



/// DECLARATIONS ET INITIALISATIONS

    // discretisation
  	int nxs[4] = {5,31,95,191}, nys[4] = {4,25,79,159};   // nxs et nys contiennent les valeurs des differents maillages
    float x_max = 1.5, y_max, h;

    // condition aux limites
    float gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    //float gs = 0.3, gw = 0.5, ge = 1., gn = 0.0, f = 0;

    // vecteurs
    float AS[MAX], AW[MAX], AP[MAX], AE[MAX], AN[MAX], b[MAX], r[MAX];

    // pour la resolution SOR
    int KMAX = 5000, K_ECH = 200;
    float EPS = 1.e-6, OMEGA = 0.6, norme2_r = 1.;
    float Q[MAX], dQ[MAX], dQ_old[MAX];

    // solution
    float T[NMAX][NMAX];

    // coordonees du maillage
    float x[NMAX], y[NMAX];

    // pour le calcul de flux
    int mx, my;
    float T1 = 293., T2 = 333., lambda = 120;
    float IH, IV, QH[NMAX], QV[NMAX];

    // pour l'etude de convergence en maillage
    int etu_conv = 1;
    float alpha, C, QV_th;

    // global
    int i, j, k, IJ;


/// POUR DIFFERENTS MAILLAGES
    for (int m = 0; m<NB_MAILLAGES ; m++){ // m pour : "maillage courant" on fait une boucle pour effectuer le programme pour plusieurs maillages

    // discretisation
    int nx = nxs[m], ny = nys[m], nn = nx * ny;           // nx  et ny  contiennent les valeurs du maillage courant
    float h = x_max / (nx + 1),y_max = (ny+1)*h;

    // Detection d'erreur de dimensions
    /*AFF*/if ((nx>NMAX)||(ny>NMAX)||(nn>MAX)){printf("ERREUR DE DIMENSIONS\nnx = %d/%d\nny = %d/%d\nnn = %d/%d\n",nx,NMAX,ny,NMAX,nn,MAX);return 0;}


/// AFFICHAGE DES PARAMETRES
    /*MAT*/fprintf(fichier_matlab,"%%%% Nouveau maillage nx = %d, ny = %d\n",nx,ny);
    /*AFF*/printf("____________________________________________\nMAILLAGE %dx%d\n--------------------------------------------\n",nx,ny);
    /*AFF*/printf("PARAMETRES\n\n");
    /*AFF*/printf("\tnn\t= %d\n\tnx\t= %d\n\tny\t= %d\n\tx_max\t= %f\n\ty_max\t= %f\n\th\t= %f\n",nn,nx,ny,x_max,y_max,h);
    /*AFF*/printf("\tgs\t= %f\n\tgw\t= %f\n\tge\t= %f\n\tgn\t= %f\n\tEps\t= %f\n",gs,gw,ge,gn,EPS);
    /*AFF*/printf("--------------------------------------------\n");

/// Mise en equation
    // initiation de la matrice pentadiagonale As, Aw, Ap, Ae, An et du second membre b
    for (i=0;i<nn;i++){// Boucle d'initialisation
        AP[i] =  4;
        AE[i] = -1;
        AW[i] = -1;
        AN[i] = -1;
        AS[i] = -1;
        b[i]  =  0;
    }
    nn = nx * ny; // sinon nn est remplace par son adresse
    for (i=0    ; i<nx ; i++  ){AS[i] = 0;b[i] += gs;}printf("AS ok\n");   // AS
    for (i=nn-nx; i<nn ; i++  ){AN[i] = 0;b[i] += gn;}printf("AN ok\n");   // AN
    for (i=nx-1 ; i<nn ; i+=nx){AE[i] = 0;b[i] += ge;}printf("AE ok\n");   // AE
    for (i=0    ; i<nn ; i+=nx){AW[i] = 0;b[i] += gw;}printf("AW ok\n");   // AW

    // impression de la matrice et du second membre
    /*AFF*/if (nn < AFF) { // si petit maillage, affichage
    /*AFF*/    printf("|AW\t|AS\t|AP\t|AN\t|AE\t|b\n");
    /*AFF*/    printf("--------------------------------------------\n");
    /*AFF*/    for (i=0;i<nn;i++){
    /*AFF*/        printf("|%2.1f\t|%2.1f\t|%2.1f\t|%2.1f\t|%2.1f\t|%2.1f\n",AW[i],AS[i],AP[i],AN[i],AE[i],b[i]);
    /*AFF*/    }
    /*AFF*/}

/// RESOLUTION SOR PAR POINT

    // initialisation
    for (i=0;i<nn+2*nx;i++){
        Q[i]  = 0.;
        dQ[i] = 0.;dQ_old[i] = 0.;
        r[i]  = 0.;
    }
    /*MAT*/fprintf(fichier_matlab,"\t%% Norm contient les iterations k sur la premiere colonne et la norme du residu sur la seconde\n\tNorm%d = [ [0 0]",m+1);
    // boucle conditionnelle
    for(k = 1;k<KMAX;k++) { // A completer
        for(i=0;i<nn;i++){
            r[i] = b[i]-(AP[i]*Q[nx+i]+AS[i]*Q[i]+AW[i]*Q[nx-1+i]+AE[i]*Q[nx+1+i]+AN[i]*Q[(2*nx)+i]);
        }
        norme2_r = norme2_vect(r, nn);
        if (norme2_r<EPS) break; // Condition d'arret

        /*AFF*/if ((k % K_ECH) == 0) {printf("k = %d  \tresidu = %e\n", k, norme2_r);}

        /*FI0*/fprintf(fichier0,"k = %d\tresidu = %e\n", k, norme2_r);//k designe ici le nombre d'iteration, norme2_r designe ici la norme du residu
        /*MAT*/fprintf(fichier_matlab,";[%d %f]", k, norme2_r);//k designe ici le nombre d'iteration, norme2_r designe ici la norme du residu

        for(i=0;i<nn;i++){
            dQ[i+nx] = OMEGA*(r[i]-AS[i]*dQ_old[i]-AW[i]*dQ_old[nx-1+i] )/AP[i];
            Q[i+nx]  += dQ[i+nx];
        }
        for(i=nx;i<nn+nx;i++){dQ_old[i] = dQ[i];} // On stocke dQ dans dQ_old
    }
    /*MAT*/fprintf(fichier_matlab," ];\n");
    /*AFF*/printf("\nconvergence en k = %d iterations, residu = %e\n\n", k, norme2_r);
    /*FI0*/fprintf(fichier0, "convergence en k = %d iterations, residu = %e\n", k, norme2_r);

    /*AFF*/if (nn < AFF) {
    /*AFF*/    printf("\n|Q\t\t|dQ\t\t|r\n");
    /*AFF*/    printf("--------------------------------------------\n");
    /*AFF*/    for (i=nx;i<nn+nx;i++){
    /*AFF*/        printf("|%f\t|%.2e\t|%.2e\n",Q[i],dQ[i],r[i-nx]);
    /*AFF*/    }
    /*AFF*/}

/// MAILLAGE ET ECRITURE DANS FICHIER_MATLAB

    /*MAT*/fprintf(fichier_matlab,"\t%% xi\n\tx%d = [",m+1);
    /*MAT*/for (i = 0; i<nx+2 ; i++){x[i]=h*i;     fprintf(fichier_matlab,"%f ",x[i]);}fprintf(fichier_matlab,"]';\n");

    /*MAT*/fprintf(fichier_matlab,"\t%% yi\n\ty%d = [",m+1);
    /*MAT*/for (i = 0; i<ny+2 ; i++){y[i]=h*i;     fprintf(fichier_matlab,"%f ",y[i]);}fprintf(fichier_matlab,"]';\n");

/// APPEL DE LA FONCTION ECRITURE

    ecriture(Q, T, x, y, nx, ny, h, gs, gw, ge, gn);

/// ECRITURE DU RESULTAT DANS FICHIER_MATLAB

    /*MAT*/fprintf(fichier_matlab,"\t%% T\n\tT%d = [",m+1);
    /*MAT*/for (j = 0; j <= ny + 1; j++) {
    /*MAT*/        for (i = 0; i <= nx + 1 ; i++) {
    /*MAT*/                fprintf(fichier_matlab,"%f ",T[i][j]);
    /*MAT*/        }
    /*MAT*/        fprintf(fichier_matlab,"; ");
    /*MAT*/}fprintf(fichier_matlab,"];\n");


/// CALCUL DU FLUX DISCRET, QV, QH

    // mx et my
    if (nx%2 == 1) {mx = (nx+1)/2.;}
    else           {etu_conv = 0;} // on ne pourra pas effectuer d'etude de convergence en maillage
    if (ny%2 == 1) {my = (ny+1)/2.;}
    else           {etu_conv = 0;} // on ne pourra pas effectuer d'etude de convergence en maillage

    QV[m] = 0.;
    QH[m] = 0.;

    if ((ny%2 == 1) && (nx%2 == 1)){ // le maillage contient donc un point central
        for (j = 0; j<ny+1 ; j++){
            QH[m] += ( T[mx+1][j]-T[mx-1][j] + T[mx+1][j+1]-T[mx-1][j+1] );
        }QH[m] = -0.25*QH[m];

        for (i = 0; i<nx+1 ; i++){
            QV[m] += ( T[i][my+1]-T[i][my-1] + T[i+1][my+1]-T[i+1][my-1] );
        }QV[m] = -0.25*QV[m];

        //printf("IH = %f, IV = %f\n", IH, IV);
        printf("QH%d = %f, QV%d = %f\n", m, QH[m], m, QV[m]);
    }
    else {printf("Le maillage n'est pas approprie au calcul de flux car nx ou ny pair\n");}

    // Stockage des 'donnees du maillage courant' dans les 'variables de sortie de la fonction matlab'
    /*MAT*/fprintf(fichier_matlab,"\n\tNorm(%d) = {Norm%d};\n",m+1,m+1);
    /*MAT*/fprintf(fichier_matlab,"\tx(%d) = {x%d};\n",m+1,m+1);
    /*MAT*/fprintf(fichier_matlab,"\ty(%d) = {y%d};\n",m+1,m+1);
    /*MAT*/fprintf(fichier_matlab,"\tT(%d) = {T%d};\n",m+1,m+1);
    /*AFF*/printf("\n\n");
    }// fin du for (int m = 0; m<3 ; m++){ // m pour : "differents maillages"


/// ETUDE DE CONVERGENCE EN MAILLAGE
    // pour cette etude, on se penche vers les trois premiers maillages

    if ((etu_conv == 1) && (NB_MAILLAGES >= 3)){// si tous les maillages (+de3) ont un point central
        h = x_max / (nxs[0] + 1); // h1

        alpha = log((QV[0]-QV[1])/(QV[1]-QV[2])) / log(2);
        C = (QV[0]-QV[1]) / ( h*h*(1-(1/pow(2,alpha))) );
        QV_th = QV[0] - C*pow(h,alpha) ;
        printf("QV_theorique = %f", QV_th);
    }
    else {printf("Etude de convergence en maillage impossible\n");}

    /*FI0*/fclose(fichier0);
    /*MAT*/fprintf(fichier_matlab,"\nend\n");
    /*MAT*/fclose(fichier_matlab);
    /*AFF*/printf("\n\nFIN DU PROGRAMME\n");

	return 0;

}
