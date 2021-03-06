
#include <NTL/fileio.h>

#include <string.h>

#include <NTL/new.h>

NTL_START_IMPL


void OpenWrite(ofstream& s, const char *name)
{
   s.open(name, ios::out);

   if (!s) {
      cerr << "open write error: " << name;
      Error("");
   }
}


void OpenRead(ifstream& s, const char *name)
{
   s.open(name, ios::in);
   if (!s) {
      cerr << "open read error: " << name;
      Error("");
   }
}

#define SBUF_SZ (1024)

NTL_THREAD_LOCAL static char sbuf[SBUF_SZ];

char *FileName(const char* stem, const char *ext)
{
   if (strlen(stem) + 1 + strlen(ext) >= SBUF_SZ) Error("bad file name");

   strcpy(sbuf, stem);
   strcat(sbuf, "-");
   strcat(sbuf, ext);

   return sbuf;
}

char *FileName(const char* stem, const char *ext, long d)
{
   if (strlen(stem) + 1 + strlen(ext) + 1 + 5 >= SBUF_SZ) Error("bad file name");

   strcpy(sbuf, stem);
   strcat(sbuf, "-");
   strcat(sbuf, ext);
   strcat(sbuf, "-");

   char dbuf[6];
   dbuf[5] = '\0';
   long i, dd;
   dd = d;
   for (i = 4; i >= 0; i--) {
      dbuf[i] = IntValToChar(dd % 10);
      dd = dd / 10;
   }

   strcat(sbuf, dbuf);

   return sbuf;
}

NTL_END_IMPL
