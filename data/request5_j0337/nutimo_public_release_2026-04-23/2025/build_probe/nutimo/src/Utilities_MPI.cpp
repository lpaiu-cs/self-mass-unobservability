// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include "Utilities_MPI.h"


/// **************** MPI utility ***************************
void MPI_Collect_and_print(const char * strtoprint, const int size_strtoprint, const int printingproc, const int * printers, const int nprinters, const MPI_Comm MPI_COMM_printers)
{
            int code;
            int k = 0;
            int mpirank = 0;
            MPI_Comm_rank(MPI_COMM_printers,&mpirank);
            char * strtoprint_remote;
            strtoprint_remote = (char*) malloc(sizeof(char) * size_strtoprint);

            if (mpirank != printingproc)
            {
                for (k = 0; k < nprinters; k++)
                {
                    if (printers[k] == mpirank)
                    {
//                         printf(">>>>>>>>>>>>>>>sending %i \n",mpirank);
                        code = MPI_Send(strtoprint, size_strtoprint, MPI_CHARACTER, printingproc, mpirank, MPI_COMM_printers);
                    }
                }
            }
            if (mpirank == printingproc)
            {
                // Collect and print information from each printer
                for (k = 0; k < nprinters; k++)
                {
                    if (printers[k] != mpirank)
                    {
//                         printf("<<<<<<<<<<<<<<<<<<receiving %i \n", printers[k]);
                        code = MPI_Recv(strtoprint_remote, size_strtoprint, MPI_CHARACTER, printers[k], printers[k], MPI_COMM_printers, MPI_STATUS_IGNORE);
                        printf("%s",strtoprint_remote);
//                         printf("<<<<<<<<<<<<<<<<<<end of receiving %i \n", printers[k]);
                    }
                    else
                    {
                        printf("%s",strtoprint);
                    }
                }

            }
        free(strtoprint_remote);
}; // End of collect and print
/// **************** end of MPI utility ***************************
