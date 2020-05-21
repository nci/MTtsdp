#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int ac, char **av) {
    // *.x   filename  nline  outfilename
    if (ac<4) 
    {
        fprintf(stderr, "Wrong. Command should be run: %s  filename  nline  outfilename\n",av[0] );
        exit(0);
    }
    char *filename = av[1];
    char *nmout = av[3];
    int nline;
    sscanf(av[2], "%d", &nline);
    FILE *fin = fopen(filename, "r");
    if (fin==NULL) {
        fprintf(stderr, "Fail in open %s\n", filename);
    }
    FILE *fout = fopen(nmout, "wb");
    if (fout==NULL) {
        fprintf(stderr, "Fail in open %s\n", nmout);
    }
    //
    float *vec = (float*) calloc(nline, sizeof(float) );
    for(int i=0; i<nline; ++i) {
        fscanf(fin, "%f", &(vec[i]) );
    }
    fwrite(vec, sizeof(float), nline, fout);
    //
    fclose(fin);
    fclose(fout);
    free(vec);
    fprintf(stdout, "Done %d\n", nline);

    return 0;
}