#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/crd.h"
struct rh1 h1;
struct rd10 d10;
struct rd30 d30;

/*-------------------------------------------------------------------------
 * Program: frd_strip
 *
 * Purpose:
 * Remove station-specific records from full rate data, and compress some
 *      other records. This is the final step to create the fullrate file
 *	to send to the data center.
 *
 * Calling sequence:
 *   frd_strip input_crd_filename output_crd_filename
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 31, 2007 - Initial version
 *   Sept 13, 2018 - For CRD format v2. (Added read and write h1 to get
 *                   and put the version number.) rlr.
 *
**-----------------------------------------------------------------------*/
int
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[256];
  FILE *str_in, *str_out;

  if ((str_in = fopen (argv[1], "r")) == NULL)
    {
      printf ("Could not open file %s\n", argv[1]);
      exit (1);
    }
  if ((str_out = fopen (argv[2], "w")) == NULL)
    {
      printf ("Could not open file %s\n", argv[2]);
      exit (1);
    }

  /* Copy and reformat data */
  while (fgets (str, 256, str_in) != NULL)
    {
      /* Read the version number from the H1 record & write it out */
      if (strncmp (str, "h1", 2) == 0 || strncmp (str, "H1", 2) == 0)
	{
	  read_h1 (str, &h1);
	  write_h1 (str_out, h1);
	}
      /* drop the station-specific records */
      else if (strncmp (str, "9", 1) == 0)
	{
	}
      /* squeeze out some spaces */
      else if (strncmp (str, "10", 2) == 0)
	{
	  read_10 (str, &d10);
	  write_10 (str_out, d10);
	}
      /* drop the extra angle records */
      else if (strncmp (str, "30", 2) == 0)
	{
	  read_30 (str, &d30);
	  if (d30.angle_origin_ind >= 0)
	    fputs (str, str_out);
	}
      /* just copy the rest */
      else
	{
	  fputs (str, str_out);
	}

    }
}
