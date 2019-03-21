/* nrerror
 * *** NUMERICAL RECIPES ERROR HANDLER ***
 *
 * This code is in the public domain.
 */

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
     printf("Numerical Recipes run-time error...\n");
     printf("%s\n",error_text);
     printf("...now exiting to system...\n");
     exit(1);
}

/* dmatrix
 * *** ALLOCATES DOUBLE PRECISION MATRICES ***
 *
 * the matrix has range m[nrl..nrh][ncl..nch]
 *
 * This code is in the public domain.
 */

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range               */
/* m[nrl..nrh][ncl..nch]                                       */
/* NR_END has been replaced with its value, 1.                 */
{
   long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += 1;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *)malloc((size_t)((nrow*ncol+1)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += 1;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* Sets the newly created matrix to zero */
   for(i=nrl;i<=nrh;i++) for(j=ncl;j<=nch;j++) m[i][j] = 0.;

   /* return pointer to array of pointers to rows */
   return m;
}

/* free_dmatrix
 * *** DE-ALLOCATES DOUBLE PRECISION MATRICES ***
 *
 * the matrix has range m[nrl..nrh][ncl..nch]
 *
 * This code is in the public domain.
 */

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free an double matrix allocated by dmatrix() */
/* replaced NR_END => 1, FREE_ARG => (char *)   */
{
   free((char *) (m[nrl]+ncl-1));
   free((char *) (m+nrl-1));
}

/* imatrix
 * *** ALLOCATES INTEGER MATRICES ***
 *
 * the matrix has range m[nrl..nrh][ncl..nch]
 *
 * This code is in the public domain.
 */

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an integer matrix with subscript range             */
/* m[nrl..nrh][ncl..nch]                                       */
/* NR_END has been replaced with its value, 1.                 */
{
   long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   int **m;

   /* allocate pointers to rows */
   m=(int **) malloc((size_t)((nrow+1)*sizeof(int*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += 1;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(int *)malloc((size_t)((nrow*ncol+1)*sizeof(int)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += 1;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* Sets the newly created matrix to zero */
   for(i=nrl;i<=nrh;i++) for(j=ncl;j<=nch;j++) m[i][j] = 0;

   /* return pointer to array of pointers to rows */
   return m;
}

/* free_imatrix
 * *** DE-ALLOCATES INTEGER MATRICES ***
 *
 * the matrix has range m[nrl..nrh][ncl..nch]
 *
 * This code is in the public domain.
 */

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an integer matrix allocated by imatrix() */
/* replaced NR_END => 1, FREE_ARG => (char *)   */
{
   free((char *) (m[nrl]+ncl-1));
   free((char *) (m+nrl-1));
}
/* End free_imatrix */

/* dvector
 * *** ALLOCATES DOUBLE PRECISION VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
/* replaced macros, as with dmatrix etc.                   */
{
   double *v;
   long i;

   v=(double *)malloc((size_t) ((nh-nl+2)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");

   /* Sets the newly created vector to zero */
   for(i=0;i<nh-nl+2;i++) v[i] = 0.;

   return(v-nl+1);
}
/* End dvector */

/* free_dvector
 * *** DE-ALLOCATES DOUBLE PRECISION VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
   free((char*) (v+nl-1));
}
/* End free_dvector */

/* lvector
 * *** ALLOCATES LONG INTEGER VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

long *lvector(long nl, long nh)
/* allocate a long vector with subscript range v[nl..nh] */
{
   long *v;
   long i;

   v=(long *)malloc((size_t) ((nh-nl+2)*sizeof(long)));
   if (!v) nrerror("allocation failure in lvector()");

   /* Sets the newly created vector to zero */
   for(i=0;i<=nh-nl;i++) v[i] = 0;

   return(v-nl+1);
}
/* End lvector */

/* free_lvector
 * *** DE-ALLOCATES LONG INTEGER VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

void free_lvector(long *v, long nl, long nh)
/* free a long vector allocated with lvector() */
{
   free((char*) (v+nl-1));
}
/* End free_lvector */

/* ivector
 * *** ALLOCATES LONG INTEGER VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

int *ivector(long nl, long nh)
/* allocate an integer vector with subscript range v[nl..nh] */
{
   int *v;

   v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
   if (!v) nrerror("allocation failure in ivector()");
   return(v-nl+1);
}
/* End ivector */

/* free_ivector
 * *** DE-ALLOCATES INTEGER VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

void free_ivector(int *v, long nl, long nh)
/* free an integer vector allocated with ivector() */
{
   free((char*) (v+nl-1));
}
/* End free_ivector */

/* NUMERICAL RECIPES NON-PUBLIC DOMAIN FUNCTION */

/* gaussjinv
 * *** INVERTS A MATRIX ***
 *
 * Takes matrix A[0..n-1][0..n-1] and inverts it.  The original matrix A is destroyed.
 *
 * Arguments:
 * > A: matrix to be inverted
 *   n: dimension of matrix
 */

void gaussjinv(double **A, int n) {
   /* This routine inverts a matrix A and stores the result back to
    * A.  The matrices are taken to be n by n.  This routine is a
    * slight modification to Numerical Recipes gaussj, see Section
    * 2.1 p.39 in NR C 2ed.  Note we've replaced float with double.
    */

   int *indxc, *indxr, *ipiv;
   int i,icol,irow,j,k,l,ll;
   double big,dum,pivinv,temp;

   /* The integer arrays ipiv, indxr, and indxc are used for    */
   /* bookkeeping on the pivoting.                              */
   indxc = ivector(0,n-1);
   indxr = ivector(0,n-1);
   ipiv = ivector(0,n-1);

   for(j=0;j<n;j++) ipiv[j]=0;
   for(i=0;i<n;i++) { /* This is the main loop over the columns */
      big=0.0;                              /* to be reduced. */
      for(j=0;j<n;j++)       /* This is the outer loop of the */
         if (ipiv[j] != 1)   /* search for a pivot element. */
            for (k=0;k<n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(A[j][k]) >= big) {
                     big=fabs(A[j][k]);
                     irow=j;
                     icol=k;
                  }
               } else if (ipiv[k] > 1)
                  nrerror("gaussj: Singular Matrix-1");
               } /* end for(k) loop */
      /* also end for(j) loop */
      ++(ipiv[icol]);

      /* We now have the pivot element, so we interchange
       * rows, if needed, to put the pivot element on the
       * diagonal.  The columns are not physically
       * interchanged, only relabeled: indxc[i], the column
       * of the ith pivot element, is the ith column that is
       * reduced, while indxr[i] is the row in which that
       * pivot element was originally located.  If indxr[i]
       * != indxc[i] there is an implied column interchange.
       * With this form of bookkeeping, the solution Bs (of
       * which we don't have any!) will end up in the correct
       * order, and the inverse matrix (i.e. that remains in
       * A) will be scrambled by columns.
       */


      if (irow != icol) {
         for (l=0;l<n;l++) {
            /* Swap A[irow][l] and A[icol][l] */
            temp = A[irow][l];
            A[irow][l] = A[icol][l];
            A[icol][l] = temp;
         }
      } /* end if */
      indxr[i]=irow; /* We are now ready to divide the pivot  */
      indxc[i]=icol; /* row by the pivot element, located at  */
                                            /* irow and icol. */
      if (A[icol][icol] == 0.0)
         nrerror("gaussj: Singular Matrix-2");
      pivinv=1.0/A[icol][icol];
      A[icol][icol]=1.0;
      for (l=0;l<n;l++) A[icol][l] *= pivinv;
      for (ll=0;ll<n;ll++) /* Now we reduce the rows, except  */
         if (ll != icol) { /* for the pivot one, of course. */
            dum=A[ll][icol];
            A[ll][icol]=0.0;
            for(l=0;l<n;l++) A[ll][l] -= A[icol][l]*dum;
         }
   }

   /* This is the end of the main loop over columns of the
    * reduction. It only remains to unscramble the solution in
    * view of the column interchanges.  We do this by
    * interchanging pairs of columns in the reverse order that
    * the permutation was built up.
    */
   for (l=n-1;l>=0;l--) {
      if (indxr[l] != indxc[l])
         for (k=0;k<n;k++) {
            temp = A[k][indxr[l]];
            A[k][indxr[l]] = A[k][indxc[l]];
            A[k][indxc[l]] = temp;
         }
   } /* And we are done. */

   /* Clean up memory */
   free_ivector(ipiv,0,n-1);
   free_ivector(indxr,0,n-1);
   free_ivector(indxc,0,n-1);
}
