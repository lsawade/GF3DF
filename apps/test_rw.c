/************************************************************

  This example shows how to read and write data to a
  dataset.  The program first writes integers to a dataset
  with dataspace dimensions of DIM0xDIM1, then closes the
  file.  Next, it reopens the file, reads back the data, and
  outputs it to the screen.

  This file is intended for use with HDF5 Library version 1.8

 ************************************************************/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#define FILE            "h5ex_d_rdwr.h5"
#define DATASET         "DS1"
#define DIM0            3
#define DIM1            3
#define DIM2            3
#define DIM3            1

int
main (void)
{
    hid_t       file, space, dset;          /* Handles */
    herr_t      status;
    hsize_t     dims[4] = {DIM0, DIM1, DIM2, DIM3};
    int         wdata[DIM0][DIM1][DIM2][DIM3],          /* Write buffer */
                rdata[DIM0][DIM1][DIM2][DIM3],          /* Read buffer */
                i, j, k, l, counter;

    /*
     * Initialize data.
     */
    counter = 0;
    for (i=0; i<DIM0; i++)
        for (j=0; j<DIM1; j++)
            for (k=0; k<DIM2; k++) {
                wdata[i][j][k][0] = counter;
                counter += 1;
            }

    printf ("%s:\n", DATASET);
    for (i=0; i<DIM0; i++) {
        printf (" [\n");
        for (j=0; j<DIM1; j++) {
            printf ("   [");
            for (k=0; k<DIM2; k++)
                printf (" %3d", wdata[i][j][k][0]);
            printf ("   ]\n");
            }
        printf ("]\n");
    }
    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple (4, dims, NULL);

    /*
     * Create the dataset.  We will use all default properties for this
     * example.
     */
    dset = H5Dcreate (file, DATASET, H5T_STD_I32LE, space, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset.
     */
    status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                wdata[0]);

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Fclose (file);


    /*
     * Now we begin the read section of this example.
     */

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen (file, DATASET, H5P_DEFAULT);

    /*
     * Read the data using the default properties.
     */
    status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                rdata[0]);

    /*
     * Output the data to the screen.
     */
    printf ("%s:\n", DATASET);
    for (i=0; i<DIM0; i++) {
        printf (" [\n");
        for (j=0; j<DIM1; j++) {
            printf ("   [");
            for (k=0; k<DIM2; k++)
                printf (" %3d", rdata[i][j][k][0]);
            printf ("   ]\n");
            }
        printf ("]\n");
    }

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Fclose (file);

    return 0;
}