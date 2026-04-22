// A C/C++ program for splitting a string using strtok() 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
  
#ifdef TEST
int main() 
{ 
    //char str[] = "Geeks-for-Geeks"; 
    char str[] = "11 11191.248834096967     0.048014200565 std 2  120.0    137  128.6   0.296  -0.520     -35.0  11.42 0"; 
    int rows= 100;
    int cols= 256;
    char *tokens = malloc(rows * cols * sizeof(char));
    int i, ntoks;
    int tokenize (char *, char *, int, int, char *);
  
    ntoks= tokenize (str, " ", 100, 256, tokens);
    for (i=0; i<ntoks; i++)
      {
        printf(">> %d [%s]\n", i, &tokens[i*cols]); 
      }
    free (tokens);
  
    return 0; 
} 
#endif

int tokenize (char *str, char *separator, int rows, int cols, char *tokens)
{
  char* token = strtok(str, separator); 
  int nt= 0;
  
  // Keep printing tokens while one of the delimiters present in str[]. 
  while (token != NULL) { 
      printf("%s\n", token); 
      strncpy (&tokens[nt*cols], token, cols); 
      tokens[nt*cols+ cols-1]= '\0';	// make sure the strin is terminated
      printf("%d [%s]\n", nt, &tokens[nt*cols]); 
      nt++;
      if (nt == rows)
	{
	  printf ("Too many tokens for array! Quitting!\n");
	  return (-1);
        }
      token = strtok(NULL, separator); 
  } 
  return (nt);
}
