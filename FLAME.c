#include <flame.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  float **data = NULL;
  char **gene_ids = NULL;
  int i, j, ki, ti, N = 0, M = 0;
  FILE *fin;
  FILE *fout;
  Flame *flame;

  int knn_vals[] = {1, 2, 4, 5, 7, 9, 11, 13};
  double membership_th_vals[] = {
    -1.0,
    0.01,
    0.1,
    0.2,
    0.3,
    0.5,
    0.7,
    0.9,
    0.95,
    0.99
  };

  int nk = sizeof(knn_vals) / sizeof(knn_vals[0]);
  int nt = sizeof(membership_th_vals) / sizeof(membership_th_vals[0]);

  if (argc < 2) {
    printf("No input file\n");
    return 0;
  }

  fin = fopen(argv[1], "r");
  if (!fin) {
    printf("Cannot open input file\n");
    return 1;
  }

  fscanf(fin, "%i %i\n", &N, &M);

  printf("Reading dataset with %i rows and %i columns\n", N, M);

  while ((i = fgetc(fin)) != '\n' && i != EOF);

  gene_ids = malloc(N * sizeof(char *));
  data = malloc(N * sizeof(float *));

  for (i = 0; i < N; i++) {
    gene_ids[i] = malloc(256 * sizeof(char));
    data[i] = malloc(M * sizeof(float));

    fscanf(fin, "%255s", gene_ids[i]);

    for (j = 0; j < M; j++) {
      fscanf(fin, "%f", &data[i][j]);
    }
  }

  fclose(fin);

  system("mkdir -p results");

  for (ki = 0; ki < nk; ki++) {
    for (ti = 0; ti < nt; ti++) {
      int knn = knn_vals[ki];
      double membership_th = membership_th_vals[ti];
      char outfile[256];

      if (membership_th < 0) {
        sprintf(outfile, "results/result_knn_%d_thresh_asterisk.txt", knn);
      } else {
        sprintf(outfile, "results/result_knn_%d_thresh_%0.2f.txt", knn, membership_th);
      }

      printf("Running knn=%d, threshold=", knn);

      if (membership_th < 0) {
        printf("*");
      } else {
        printf("%g", membership_th);
      }

      printf(" ... ");

      flame = Flame_New();

      Flame_SetDataMatrix(flame, data, N, M, 0);
      Flame_DefineSupports(flame, knn, -2.0);
      Flame_LocalApproximation(flame, 500, 1e-6);
      Flame_MakeClusters(flame, membership_th);

      fout = fopen(outfile, "w");
      if (!fout) {
        printf("Cannot open output file\n");
        Flame_Clear(flame);
        free(flame);
        continue;
      }

      fprintf(fout, "Dataset: %i rows and %i columns\n", N, M);
      fprintf(fout, "knn = %i\n", knn);
      fprintf(fout, "support_density_threshold = -2.0\n");

      if (membership_th < 0) {
        fprintf(fout, "membership_threshold = *\n");
      } else {
        fprintf(fout, "membership_threshold = %g\n", membership_th);
      }

      fprintf(fout, "Found %i cluster supporting objects\n", flame->cso_count);
      fprintf(fout, "Final clusters including outlier group: %i\n\n", flame->count);

      for (i = 0; i < flame->count; i++) {
        if (i == flame->count - 1) {
          fprintf(fout, "\nCluster outliers, with %6i members:\n",
                  flame->clusters[i].size);
        } else {
          fprintf(fout, "\nCluster %3i, with %6i members:\n",
                  i + 1, flame->clusters[i].size);
        }

        for (j = 0; j < flame->clusters[i].size; j++) {
          int gene_index = flame->clusters[i].array[j];

          if (j) {
            fprintf(fout, ",");
            if (j % 10 == 0) fprintf(fout, "\n");
          }

          fprintf(fout, "%s", gene_ids[gene_index]);
        }

        fprintf(fout, "\n");
      }

      fclose(fout);

      Flame_Clear(flame);
      free(flame);

      printf("done\n");
    }
  }

  for (i = 0; i < N; i++) {
    free(data[i]);
    free(gene_ids[i]);
  }

  free(data);
  free(gene_ids);

  printf("All done. Results are in the results/ directory.\n");

  return 0;
}