#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
	return findMatrixValue (x*nx, y*ny, z*nz);
}

double normaVet (double vet[3])
{
	return sqrt(vet[0]*vet[0] + vet[1]*vet[1] + vet[2]*vet[2]);
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
	ajustaTamanhoVet(dir, -sqrt(3.0));

	ponto[0] += 0.5;
	ponto[1] += 0.5;
	ponto[2] += 0.5;
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
	
	fprintf(arq, "P2\n# %s\n%d %d\n%d\n", nome, lin, col, normalize_char(maior, maior));
	printf("P2\n# %s\n%d %d\n%d\n", nome, lin, col, normalize_char(maior, maior));
	for (i=0; i<lin; i++)
	{
		for (j=0; j<col; j++)
		{
			fprintf(arq, "%d, ", normalize_char(mat[i][j], maior));
			//printf("%d ", normalize_char(mat[i][j], maior));
		}
		fprintf(arq, "\n");
		//printf("\n");
	}

}

double DoubleSimpson (double a[3], double b[3], double* val)
{
	double Sab, Sac, Scb, Eab_16, h = norma(b-a), c=midpoint(a,b), 
		fa=findMatrixValueNormal(a), fb=findMatrixValueNormal(b), fc = findMatrixValueNormal(c), fm = fc;

	Sab = (h/6)*(fa+4*fm+fb);

	Sac = (h/12)*(fa+4*findMatrixValueNormal(midpoint(a,c))+fc);
	Scb = (h/12)*(fc+4*findMatrixValueNormal(midpoint(c,b))+fb);

	Eab_16 = (Sab - Sac - Scb)/15;

	if (Eab_16<0)
		Eab_16 = -Eab_16;

	*val = Sac + Scb + Eab_16;

	return Eab_16;
}

double AdaptiveSimpson (double a, double b, double (*f) (double x), double tol)
{
	double diff, resp, c=(a+b)/2;

	diff = DoubleSimpson (a, b, f, &resp);

	if (diff<tol)
		return resp;

	return AdaptiveSimpson(a, c, f, tol/2) + AdaptiveSimpson(c, b, f, tol/2); 
}

//recebe a direcao, retorna todo o resto
void preparaParaInicio (double dir[3], double plano[4], double direita[3], double cima[3], double pontoinicial[3], double tampasso)
{
	int i;
	double meio[3], sqrt3 = sqrt(3.0);
	acharPontoNaEsfera (dir, meio);
	acharPlanoTangenteNoPonto (meio, plano);

	//Um vetor qualquer no plano, que passará a ser o nosso conceito de cima
	cima[0] = -meio[0];
	cima[1] = -meio[1];
	cima[2] = plano[3]/plano[2] - meio[2];

	prodVetorial (dir, cima, direita);

	ajustaTamanhoVet(cima, tampasso);
	ajustaTamanhoVet(direita, tampasso);
	ajustaTamanhoVet(dir, 1);

	for (i=0; i<3; i++)
		pontoinicial[i] = meio[i] - sqrt3/2*cima[i] -sqrt3/2*direita[i];
}

unsigned char** geraImagemComIntegralVersao1 (double dir[3], double h, int tamfinal)
{
	double plano[4], direita[3], cima[3], pontoinicial[3], ponto[3];
	double passo = 1/tamfinal;
	int i, j;
	unsigned char** img;

	img = (unsigned char**) malloc (sizeof(unsigned char*)*tamfinal);
	for (i=0; i<tamfinal; i++)
		img[i] = (unsigned char*) malloc (sizeof(unsigned char)*tamfinal);

	preparaParaInicio(dir, plano, direita, cima, pontoinicial, passo);

	ponto[0] = pontoinicial[0];
	ponto[1] = pontoinicial[1];
	ponto[2] = pontoinicial[2];
	
	for (i=0; i<tamfinal; i++)
	{
		for (j=0; j<tamfinal; j++)
		{
			img[i][j] = funcintegra(ponto, dir, h);

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
	unsigned char **resp, maior;

	ini = (unsigned char*) malloc (nx*ny*nz);
	arq = fopen("head-8bit.raw", "rb");

	fread(ini, sizeof(unsigned char), nx*ny*nz, arq);

	resp = gera_imagem_com_media (&maior);

	gera_PGM(resp, nx, nz, maior);

	//print_matriz(resp, NX, NZ);

	fclose(arq);
	free(ini);
	for (i=0; i<NX; i++)
		free(resp[i]);
	free(resp);
	return 0;
}