#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#if defined(_WIN32) && defined(__INTEL_COMPILER)
#  include "dirent_windows.h"
#  include <io.h>
#  define F_OK 0
#  define access _access
#else
#  include <unistd.h>
#  include <dirent.h>
#endif


void file_info(const char*filename,int*mode,int*exist,int*time){
  int k;
  struct stat buf;
  k=stat(filename,&buf);
  if(k != 0) {
    *mode=0;
    *exist=0;
    *time=0;
  }else{
    *mode=buf.st_mode;
    if(*mode == 0) *exist=0; else *exist=1;
    *time=buf.st_mtime;
  }
}

void dir_exists_c(const char*dirname, int*exists)
{
    DIR* dir = opendir(dirname);

    if (dir) {
        /* Directory exists. */
      closedir(dir);
      *exists = 1;
      printf("%d", *exists);
      fflush(stdout);
    } else if (ENOENT == errno) {
    /* Directory does not exist. */
      *exists = 0;
    } else {
    /* opendir() failed for some other reason. */
      *exists = 0;
    }
}