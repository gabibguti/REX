#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TOL 0.1
#define PASSOERRO 0.01
#define RESOLUCAO 1000
#define NUMPASSOS 1000

int niteracao = 0;

int nx = 256;
int ny = 256;
int nz = 99;
unsigned char *ini;

//testar
double findMatrixValue (double x, double y, double z)
{
	int underx = (int) x, undery = (int) y, underz = (int) z,
		upperx = underx + 1, uppery = undery + 1, upperz = underz+1;
	
	double remainx = x - floor(x), remainy = y - floor(y), remainz = z - floor(z);

	double corner000 = (double) ini[underx  + undery*nx + underz*nx*ny],
		   corner001 = (double) ini[underx  + undery*nx + upperz*nx*ny],
		   corner010 = (double) ini[underx  + uppery*nx + underz*nx*ny],
		   corner011 = (double) ini[underx  + uppery*nx + upperz*nx*ny],
		   corner100 = (double) ini[upperx  + undery*nx + underz*nx*ny],
		   corner101 = (double) ini[upperx  + undery*nx + upperz*nx*ny],
		   corner110 = (double) ini[upperx  + uppery*nx + underz*nx*ny],
		   corner111 = (double) ini[upperx  + uppery*nx + upperz*nx*ny];

	double corner00 = corner000*remainz + corner001*(1-remainz),
		   corner01 = corner010*remainz + corner011*(1-remainz),
		   corner10 = corner100*remainz + corner101*(1-remainz),
		   corner11 = corner110*remainz + corner111*(1-remainz);

	double corner0 = corner00*remainy + corner01*(1-remainy),
		   corner1 = corner10*remainy + corner11*(1-remainy);

	double resp = corner0*remainx + corner1*(1-remainx);

	return resp;
}

//testar
double findMatrixValueNormal (double vet[3])
{
	double x = vet[0], y = vet[1], z = vet[2];
	if (x<0||y<0||z<0||x>1||y>1||z>1)
		return 0;
	return findMatrixValue (x*nx-1, y*ny-1, z*nz-1);
}

double normaVet (double vet[3])
{
	return sqrt(vet[0]*vet[0] + vet[1]*vet[1] + vet[2]*vet[2]);
}

double normaDifVet (double vet1[3], double vet2[3])
{
	return sqrt((vet1[0]-vet2[0])*(vet1[0]-vet2[0]) + (vet1[1]-vet2[1])*(vet1[1]-vet2[1]) + (vet1[2]-vet2[2])*(vet1[2]-vet2[2]));
}

void ajustaTamanhoVet (double vet[3], double newsize)
{
	double tam = normaVet(vet);
	vet[0] = vet[0]/tam*newsize;
	vet[1] = vet[1]/tam*newsize;
	vet[2] = vet[2]/tam*newsize;
}

void acharPontoNaEsfera (double dir[3], double ponto[3])
{
	ajustaTamanhoVet(dir, -sqrt(3.0)/2.);

	ponto[0] = dir[0] + 0.5;
	ponto[1] = dir[1] + 0.5;
	ponto[2] = dir[2] + 0.5;

	ajustaTamanhoVet (dir, -sqrt(3.0)); //tamanho 1, com a direção inicial
}

void acharPlanoTangenteNoPonto (double ponto[3], double plano[4])
{
	plano[0] = 0.5 - ponto[0];
	plano[1] = 0.5 - ponto[1];
	plano[2] = 0.5 - ponto[2];

	plano[3] = 0.5*(ponto[0]+ponto[1]+ponto[2]) - (ponto[0]*ponto[0]+ponto[1]*ponto[1]+ponto[2]*ponto[2]);
}

void prodVetorial (double v1[3], double v2[3], double resp[3])
{
	resp[0] = v1[1]*v2[2]-v1[2]*v2[1];
	resp[1] = v1[2]*v2[0]-v1[0]*v2[2];
	resp[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

unsigned char normalize_char (unsigned char num, unsigned char max)
{
	double temp = (double) num;
	temp = temp*255./(double)max;
	return (unsigned char) temp;
}

double cvtdouble (unsigned char x)
{
	return (double)x;
}

unsigned char cvtchar (double x)
{
	return (unsigned char)x;
}

void print_matriz (unsigned char** mat, int lin, int col)
{
	int i, j;
	for (i=0; i<lin; i++)
	{
		for (j=0; j<col; j++)
			printf("%02x ", mat[i][j]);
		printf("\n");
	}
	printf("\n");
}

unsigned char** gera_imagem_com_media (unsigned char *maior)
{
	int i, j, k;
	unsigned char **resp;
	double soma;

	*maior=0;

	resp = (unsigned char**) malloc (sizeof(unsigned char*)*nx);
	for (i=0; i<nx; i++)
		resp[i] = (unsigned char*) malloc (sizeof(unsigned char)*nz);

	for (i=0; i<nx; i++)
	{
		for (k=0; k<nz; k++)
		{
			soma=0;
			for (j=0; j<ny; j++)
			{
				soma += cvtdouble(ini[i  + j*nx + k*nx*ny]);
			}
			if ((resp[i][k] = cvtchar(soma/ny))>*maior)
				*maior = resp[i][k];
		}
	}



	return resp;
}

void gera_PGM (unsigned char** mat, int lin, int col, unsigned char maior)
{
	int i, j;
	char nome[] = "caveira.csv";
	FILE* arq = fopen (nome, "w");
	
	//fprintf(arq, "P2\n# %s\n%d %d\n%d\n", nome, lin, col, normalize_char(maior, maior));
	printf("P2\n# %s\n%d %d\n%d\n", nome, lin, col, normalize_char(maior, maior));
	for (i=0; i<lin; i++)
	{
		for (j=0; j<col; j++)
		{
			fprintf(arq, "%d, ", mat[i][j]);
			//printf("%d ", normalize_char(mat[i][j], maior));
		}
		fprintf(arq, "\n");
		//printf("\n");
	}

}

void midpoint (double v1[3], double v2[3], double v3[3])
{
	int i;
	for (i=0; i<3; i++)
		v3[i] = (v1[i]+v2[i])/2.;

}

double DoubleSimpson (double a[3], double b[3], double c[3], double* val)
{
	int i;
	double Sab, Sac, Scb, Eab_16, ab[3], h, m_ac[3], m_cb[3],
		fa=findMatrixValueNormal(a), fb=findMatrixValueNormal(b), fc = findMatrixValueNormal(c),
		fm_ab = fc, fm_ac, fm_cb;

	midpoint(a,c, m_ac);
	midpoint(c,b, m_cb);
	fm_ac = findMatrixValueNormal(m_ac);
	fm_cb = findMatrixValueNormal(m_cb);

	a[0];a[1];a[2];b[0];b[1];b[2];c[0];c[1];c[2];

	for (i=0; i<3; i++)
		ab[i] = b[i]-a[i];
	h = normaVet(ab);


	Sab = (h/6)*(fa+4*fm_ab+fb);

	Sac = (h/12)*(fa+4*fm_ac+fc);
	Scb = (h/12)*(fc+4*fm_cb+fb);

	Eab_16 = (Sab - Sac - Scb)/15;

	if (Eab_16<0)
		Eab_16 = -Eab_16;

	*val = Sac + Scb + Eab_16;

	return Eab_16;
}

double AdaptiveSimpson (double a[3], double b[3], double tol)
{
	double diff, resp, c[3];
	midpoint(a,b,c);

	diff = DoubleSimpson (a, b, c, &resp);

	if (diff<tol || niteracao>7)
		return resp;
	
	niteracao++;
	return AdaptiveSimpson(a, c, tol/2) + AdaptiveSimpson(c, b, tol/2); 
}

void vetQualquerNoPlano (double plano[4], double meio[3], double vet[3])
{
	if (plano[2]<=-1E-10 || plano[2]>=1E-10) //z =/= 0
	{
		vet[0] = 0					- meio[0];
		vet[1] = 0					- meio[1];
		vet[2] = plano[3]/plano[2]	- meio[2];
	}
	else if (plano[1]<=-1E-10 || plano[1]>=1E-10)// y =/= 0
	{
		vet[0] = 0					- meio[0];
		vet[1] = plano[3]/plano[1]	- meio[1];
		vet[2] = 0					- meio[2];
	}
	else if (plano[0]<=-1E-10 || plano[0]>=1E-10)
	{
		vet[0] = plano[3]/plano[0]	- meio[0];
		vet[1] = 0					- meio[1];
		vet[2] = 0					- meio[2];
	}
	else
		printf("Erro: Vetor direcao problematico\n");

}

//recebe a direcao, retorna todo o resto
void preparaParaInicio (double dir[3], double plano[4], double direita[3], double cima[3], double pontoinicial[3], double tampasso)
{
	int i;
	double meio[3], sqrt3 = sqrt(3.0);
	acharPontoNaEsfera (dir, meio);
	acharPlanoTangenteNoPonto (meio, plano);

	//Um vetor qualquer paralelo ao plano, que passara a ser o nosso conceito de cima
	vetQualquerNoPlano(plano, meio, cima);

	meio[0];	meio[1];	meio[2];

	prodVetorial (dir, cima, direita);

	ajustaTamanhoVet(cima, 1);
	ajustaTamanhoVet(direita, 1);

	for (i=0; i<3; i++)
		pontoinicial[i] = meio[i] - (sqrt3/2.)*cima[i] - (sqrt3/2.)*direita[i];

	ajustaTamanhoVet(cima, tampasso);
	ajustaTamanhoVet(direita, tampasso);
}

unsigned char valMedio (double ini[3], double dir[3])
{
	int i;
	double a[3], b[3], ab[3], sizeab, tam, integral, valmat, sqrt3 = sqrt(3.0), sizea, sizeb;

	for (i=0; i<3; i++)
	{
		a[i] = ini[i];
		b[i] = a[i] + dir[i];
	}
	/*
	valmat = findMatrixValueNormal(a);
	sizea = normaDifVet(a, ini);
	while (valmat==0 && sizea<=sqrt3) //Encontrar o primeiro valor valido da matriz
	{
		for (i=0; i<3; i++)
			a[i] = a[i] + dir[i]*PASSOERRO;
		valmat = findMatrixValueNormal(a);
		sizea = normaDifVet(a, ini);
	}
	
	if (valmat==0)
		return 0;

	valmat = findMatrixValueNormal(b);
	sizeb = normaDifVet(b, ini);
	while(valmat==0 && sizeb<=sqrt3) //Encontrar o ultimo valor valido da matriz
	{
		for (i=0; i<3; i++)
			b[i] = b[i] - dir[i]*PASSOERRO;
		valmat = findMatrixValueNormal(b);
		sizeb = normaDifVet(b, ini);
	}
	
	if (valmat==0)
		return 0;
	
	for (i=0; i<3; i++)
		ab[i] = b[i] - a[i]; 
	*/
	sizeab = normaDifVet(a, b);

	niteracao = 0;
	a[0]; a[1]; a[2]; b[0]; b[1]; b[2];
	integral = AdaptiveSimpson (a, b, TOL);

	return (unsigned char) (integral/sizeab); //(valor da integral)/
}

unsigned char** geraImagemComIntegralVersao1 (double dir[3], int tamfinal)
{
	double plano[4], direita[3], cima[3], pontoinicial[3], ponto[3];
	double passo = sqrt(3.)/tamfinal;
	int i, j;
	unsigned char** img;

	img = (unsigned char**) malloc (sizeof(unsigned char*)*tamfinal);
	for (i=0; i<tamfinal; i++)
		img[i] = (unsigned char*) malloc (sizeof(unsigned char)*tamfinal);

	preparaParaInicio(dir, plano, direita, cima, pontoinicial, passo);

	ponto[0] = pontoinicial[0];
	ponto[1] = pontoinicial[1];
	ponto[2] = pontoinicial[2];
	
	cima[0];cima[1];cima[2];
	direita[0];direita[1];direita[2];

	for (i=0; i<tamfinal; i++)
	{
		for (j=0; j<tamfinal; j++)
		{
			img[i][j] = valMedio (ponto, dir);

			if (img[i][j]!=0)
				printf("val\n");

			//ir para a direita
			ponto[0] += direita[0];
			ponto[1] += direita[1];
			ponto[2] += direita[2];
		}

		//resetar a direita e ir para cima
		ponto[0] = pontoinicial[0]+cima[0]*(double)i;
		ponto[1] = pontoinicial[1]+cima[1]*(double)i;
		ponto[2] = pontoinicial[2]+cima[2]*(double)i;
	}

	return img;
}

unsigned char Media (double ini[3], double dir[3], int npassos)
{
	int i;
	double soma=0, ponto[3];
	ajustaTamanhoVet(dir, sqrt(3.0)/npassos);

	ponto[0]=ini[0];
	ponto[1]=ini[1];
	ponto[2]=ini[2];

	for (i=0; i<npassos; i++)
	{
		soma += findMatrixValueNormal(ponto);
		
		ponto[0]+=dir[0];
		ponto[1]+=dir[1];
		ponto[2]+=dir[2];

	}

	return (unsigned char) (soma/ ((double) npassos));
}

unsigned char** geraImagemComMediaVersao2 (double dir[3], int tamfinal)
{
	double plano[4], direita[3], cima[3], pontoinicial[3], ponto[3];
	double passo = sqrt(3.)/tamfinal;
	int i, j;
	unsigned char** img;

	img = (unsigned char**) malloc (sizeof(unsigned char*)*tamfinal);
	for (i=0; i<tamfinal; i++)
		img[i] = (unsigned char*) malloc (sizeof(unsigned char)*tamfinal);

	preparaParaInicio(dir, plano, direita, cima, pontoinicial, passo);

	ponto[0] = pontoinicial[0];
	ponto[1] = pontoinicial[1];
	ponto[2] = pontoinicial[2];
	
	cima[0];cima[1];cima[2];
	direita[0];direita[1];direita[2];

	for (i=0; i<tamfinal; i++)
	{
		for (j=0; j<tamfinal; j++)
		{
			img[i][j] = Media (ponto, dir, NUMPASSOS);

			//ir para a direita
			ponto[0] += direita[0];
			ponto[1] += direita[1];
			ponto[2] += direita[2];
		}

		//resetar a direita e ir para cima
		ponto[0] = pontoinicial[0]+cima[0]*(double)i;
		ponto[1] = pontoinicial[1]+cima[1]*(double)i;
		ponto[2] = pontoinicial[2]+cima[2]*(double)i;
	}

	return img;
}
int main (void)
{
	int i, j, k;
	FILE* arq;
	unsigned char **resp, maior=255;
	double dir[3]={0,1,0};

	ini = (unsigned char*) malloc (nx*ny*nz);
	arq = fopen("head-8bit.raw", "rb");

	fread(ini, sizeof(unsigned char), nx*ny*nz, arq);

	resp = geraImagemComMediaVersao2 (dir, RESOLUCAO);

	gera_PGM(resp, RESOLUCAO, RESOLUCAO, maior);

	//print_matriz(resp, NX, NZ);

	fclose(arq);
	free(ini);
	for (i=0; i<RESOLUCAO; i++)
		free(resp[i]);
	free(resp);
	return 0;
}