/*
  John M. Gaspar (jsh58@unh.edu)
  June 2013 (updated 1/14, 3/14)

  This program can both filter and denoise reads
    produced by 454 pyrosequencing.

  For a complete description of the usage and
    parameters, please see the accompanying README.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FlowClus.h"

// global variables
static char* line;
static Primer* primo;
void freeMemory(void);

/* void usage(void)
 * Prints usage to stderr.
 */
void usage(void) {
  fprintf(stderr, "\n*** For a complete description of the parameters, see the README ***\n");
  fprintf(stderr, "\nUsage: ./FlowClus %s master.csv [optional parameters]\n\n", MASTERFILE);
  fprintf(stderr, "Required parameters:\n");
  fprintf(stderr, "  %s   Input master file with primer and mid tag sequences\n", MASTERFILE);
  fprintf(stderr, "Optional parameters:\n");
  fprintf(stderr, "  %s  Option to print status updates while running\n", STATUSOPT);
  fprintf(stderr, "  %s   Option to perform filtering only\n", CLEANOPT);
  fprintf(stderr, "  %s   Option to perform denoising only\n", DENOPT);
  fprintf(stderr, "  %s  Option to perform filtering and denoising\n", BOTHOPT);
  fprintf(stderr, "  %s   Input sff.txt file (req'd if filtering)\n", SFFFILE);
  fprintf(stderr, "  %s   File extension for filtered flowgrams\n", FLOWEXT);
  fprintf(stderr, "  %s   Output fasta file after filtering\n", OUTFILE);
  fprintf(stderr, "  %s   Output fasta file after denoising\n", DENFASTA);
  fprintf(stderr, "  %s   Option to produce \"QIIME-style\" output fasta file(s)\n", NOMIDOPT);
  fprintf(stderr, "  %s   Output file for counts of truncated and eliminated reads\n", ERRFILE);
  fprintf(stderr, "  %s  Output file for detailed filtering information for each read\n", FILFILE);
  fprintf(stderr, "  %s   Output file for denoising \"misses\"\n", MISSFILE);
  fprintf(stderr, "  %s   Option to produce consensus flowgram and mapping files after denoising\n", DENPOPT);
  fprintf(stderr, "  %s  File extension for denoised flowgrams\n", DENFEXT);
  fprintf(stderr, "  %s  File extension for mapping files\n", DENMEXT);
  fprintf(stderr, "  %s  Option to produce output fasta files for de novo chimera checking\n", CHIMOPT);
  fprintf(stderr, "  %s  File extension for output fasta files for UCHIME\n", UCHEXT);
  fprintf(stderr, "  %s  File extension for output fasta files for Perseus\n", PEREXT);
  fprintf(stderr, "  %s  File extension for output mapping file for chimera checking\n", UMAPEXT);
  fprintf(stderr, "  %s  Input file containing distances for each flow value\n", SDFILE);
  fprintf(stderr, "  %s  Mismatches to mid tag sequence to allow\n", MIDMIS);
  fprintf(stderr, "  %s  Mismatches to primer sequence to allow\n", PRIMMIS);
  fprintf(stderr, "  %s   Minimum sequence length\n", MINSLEN);
  fprintf(stderr, "  %s   Maximum sequence length for elimination\n", MAXSLEN);
  fprintf(stderr, "  %s   Maximum sequence length for truncation\n", MAXTRLEN);
  fprintf(stderr, "  %s   Maximum number of ambiguous bases to allow\n", MAXAMBIG);
  fprintf(stderr, "  %s   Maximum number of ambiguous bases to allow before truncation\n", OKAMBIG);
  fprintf(stderr, "  %s   Maximum homopolymer length to allow\n", MAXHOMO);
  fprintf(stderr, "  %s   Maximum homopolymer length to allow before truncation\n", OKHOMO);
  fprintf(stderr, "  %s   Option to remove the reverse primer from reads\n", REVMOPT);
  fprintf(stderr, "  %s  Option to require the reverse primer in reads\n", REVQOPT);
  fprintf(stderr, "  %s  Mismatches to reverse primer sequence to allow\n", REVMIS);
  fprintf(stderr, "  %s   Minimum average quality score\n", AVGQUAL);
  fprintf(stderr, "  %s  Length of sliding window of quality scores\n", WINDOWLEN);
  fprintf(stderr, "  %s  Minimum average quality score for sliding window\n", WINDOWAVG);
  fprintf(stderr, "  %s  Option to eliminate a read with a low quality window\n", WINDOWOPT);
  fprintf(stderr, "  %s   Absolute maximum flow value to analyze\n", MAXFLOW);
  fprintf(stderr, "  %s  Minimum flowgram length\n", MINFLEN);
  fprintf(stderr, "  %s  Maximum flowgram length\n", MAXFLEN);
  fprintf(stderr, "  %s   Noisy interval flow value minimum\n", MININT);
  fprintf(stderr, "  %s   Noisy interval flow value maximum\n", MAXINT);
  fprintf(stderr, "  %s   Maximum flow value (for truncation)\n", MAXNVAL);
  fprintf(stderr, "  %s   Minimum flow value for 4 consecutive flows\n", NOFLOW);
  fprintf(stderr, "  %s   Constant value for denoising\n", CINTER);
  fprintf(stderr, "  %s   Number of distances for denoising\n", ZINTER);
  fprintf(stderr, "  %s  Option to denoise using a trie\n", TRIEOPT);
  fprintf(stderr, "\n");

  freeMemory();
  exit(-1);
}

/* int error(char*)
 * Prints the given error message to stderr.
 */
int error(char* val, int err) {
  fprintf(stderr, "Error!  %s", val);
  char* msg;
  if (err == ERRMEM)
    msg = MERRMEM;
  else if (err == ERRCLOSE)
    msg = MERRCLOSE;
  else if (err == ERRPARAM)
    msg = MERRPARAM;
  else if (err == ERRFLOAT)
    msg = MERRFLOAT;
  else if (err == ERRINT)
    msg = MERRINT;
  else if (err == ERROPENR)
    msg = MERROPENR;
  else if (err == ERROPENW)
    msg = MERROPENW;
  else if (err == ERRINVAL)
    msg = MERRINVAL;
  else if (err == ERREXIST)
    msg = MERREXIST;
  else if (err == ERRLOAD)
    msg = MERRLOAD;
  else if (err == ERRNEED)
    msg = MERRNEED;
  else if (err == ERRSLEN)
    msg = MERRSLEN;
  else if (err == ERRWIND)
    msg = MERRWIND;
  else if (err == ERRINTER)
    msg = MERRINTER;
  else if (err == ERRFLEN)
    msg = MERRFLEN;
  else if (err == ERRMAXFL)
    msg = MERRMAXFL;
  else if (err == ERRSDVAL)
    msg = MERRSDVAL;
  else if (err == ERRDEN)
    msg = MERRDEN;
  else if (err == ERRDENVAL)
    msg = MERRDENVAL;
  else if (err == ERRPRIM)
    msg = MERRPRIM;
  else if (err == ERRPREP)
    msg = MERRPREP;
  else if (err == ERRMREP)
    msg = MERRMREP;
  else if (err == ERRNFLOWS)
    msg = MERRNFLOWS;
  else if (err == ERRHEAD)
    msg = MERRHEAD;
  else if (err == ERRORDER)
    msg = MERRORDER;
  else if (err == ERRMID)
    msg = MERRMID;
  else if (err == ERRREV)
    msg = MERRREV;
  else if (err == ERRNOREV)
    msg = MERRNOREV;
  else if (err == ERRMISM)
    msg = MERRMISM;
  else if (err == ERRMAXF)
    msg = MERRMAXF;
  else if (err == ERRLEN)
    msg = MERRLEN;
  else
    msg = UNKNOWN;

  fprintf(stderr, "%s\n", msg);
  if (err != ERREXIST && err != ERRLEN)
    freeMemory();
  return -1;
}

/* void freeTrie(Node*)
 * Recursively frees the trie.
 */
void freeTrie(Node* n) {
  if (n == NULL)
    return;
  Read* r = n->first;
  Read* tempr;
  while (r != NULL) {
    free(r->header);
    tempr = r;
    r = r->next;
    free(tempr);
  }
  if (n->flow != NULL)
    free(n->flow);
  freeTrie(n->child);
  freeTrie(n->next);
  free(n);
}

/* void freeClus(Cluster*)
 * Frees clusters and reads.
 */
void freeClus(Cluster* c) {
  Cluster* tempc;
  while (c != NULL) {
    Read* r = c->first;
    Read* tempr;
    while (r != NULL) {
      free(r->header);
      tempr = r;
      r = r->next;
      free(tempr);
    }
    free(c->flows);
    free(c->weight);
    tempc = c;
    c = c->next;
    free(tempc);
  }
}

/* void freeMemory(void)
 * Frees the stored pointers.
 */
void freeMemory(void) {
  Primer* p = primo;
  Primer* temp;
  while (p != NULL) {
    free(p->name);
    free(p->seq);
    if (p->rev != NULL)
      free(p->rev);
    free(p->dummy);
    Midtag* m = p->first;
    Midtag* tempm;
    while (m != NULL) {
      free(m->name);
      free(m->seq);
      tempm = m->next;
      free(m);
      m = tempm;
    }

    freeClus(p->head);
    freeTrie(p->root);

    if ((p->out != NULL && fclose(p->out)) ||
        (p->den != NULL && fclose(p->den)) ||
        (p->map != NULL && fclose(p->map)) ||
        (p->per != NULL && fclose(p->per)) ||
        (p->cmap != NULL && fclose(p->cmap)))
      exit(error("", ERRCLOSE));

    temp = p->next;
    free(p);
    p = temp;
  }

  free(line);
}

/* void* memalloc(int)
 * Returns pointer to allocated memory block.
 */
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* float getFloat(char*)
 * Converts the given char* to a float.
 */
float getFloat(char* in) {
  char** endptr = NULL;
  float ans = strtof(in, endptr);
  if (endptr != '\0')
    exit(error(in, ERRFLOAT));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char** endptr = NULL;
  int ans = (int) strtol(in, endptr, 10);
  if (endptr != '\0')
    exit(error(in, ERRINT));
  return ans;
}

/* int clusSize(Primer*)
 * Returns the number of clusters for this primer.
 */
int clusSize(Primer* p) {
  int clus = 0;
  for (Cluster* c = p->head; c != NULL; c = c->next)
    clus++;
  return clus;
}

/* void mergeSort()
 * Sorts the clusters based on size.
 */
void mergeSort(Primer* p, Cluster* dummy) {
  dummy->next = p->head;
  int len = clusSize(p);

  Cluster* last;
  for (int i = 1; i < len; i *= 2) {
    last = dummy;
    Cluster* c = last->next;
    Cluster* d = c;

    while (c != NULL) {
      int csize = 0;
      int dsize = i;

      // move d down the list
      for (int k = i; k; k--) {
        d = d->next;
        csize++;
        if (d == NULL) {
          dsize = 0;
          break;
        }
      }

      // build new list
      while (csize || dsize) {
        if (!dsize || (csize && c->weight[0] >= d->weight[0])) {
          last->next = c;
          last = c;
          csize--;
          c = c->next;
        } else {
          last->next = d;
          last = d;
          dsize--;
          d = d->next;
          if (d == NULL)
            dsize = 0;
        }
      }

      // reset c
      c = d;
    }

    last->next = NULL;
  }

  p->head = dummy->next;
  p->tail = last;
}

/* int clusSize2(Cluster*)
 * Returns the number of reads in this cluster.
 */
int clusSize2(Cluster* c) {
  int read = 0;
  for (Read* r = c->first; r != NULL; r = r->next)
    read++;
  return read;
}

/* void mergeSort2()
 * Sorts the reads of a cluster based on size.
 */
void mergeSort2(Cluster* c, Read* dummy) {
  dummy->next = c->first;
  int len = clusSize2(c);

  Read* last;
  for (int i = 1; i < len; i *= 2) {
    last = dummy;
    Read* c = last->next;
    Read* d = c;

    while (c != NULL) {
      int csize = 0;
      int dsize = i;

      // move d down the list
      for (int k = i; k; k--) {
        d = d->next;
        csize++;
        if (d == NULL) {
          dsize = 0;
          break;
        }
      }

      // build new list
      while (csize || dsize) {
        if (!dsize || (csize && c->length <= d->length)) {
          last->next = c;
          last = c;
          csize--;
          c = c->next;
        } else {
          last->next = d;
          last = d;
          dsize--;
          d = d->next;
          if (d == NULL)
            dsize = 0;
        }
      }

      // reset c
      c = d;
    }

    last->next = NULL;
  }

  c->first = dummy->next;
}

/* void printMiss(FILE*)
 * Prints the "miss" array to the given file.
 */
void printMiss(FILE* out, int** miss, int max) {

  // print array
  for (int i = 0; i < max; i++) {
    fprintf(out, "%d", miss[i][0]);
    for (int j = 1; j < max; j++)
      fprintf(out, ",%d", miss[i][j]);
    fprintf(out, "\n");
  }

  // free array
  for (int i = 0; i < max; i++)
    free(miss[i]);
  free(miss);

  if (fclose(out))
    exit(error("", ERRCLOSE));
}

/* void printPers()
 * Prints to the fasta file that can be further
 *   analyzed for PCR chimeras by UCHIME or Perseus.
 */
void printPers(Primer* p, Cluster* c, char* seq,
    int count, int clus, int perOpt) {

  // print to fasta file
  if (perOpt)
    fprintf(p->per, ">%s%s%d%s%d\n%s\n", c->lon->header,
      PER, clus, PER, count, seq);
  else
    fprintf(p->per, ">%s /ab=%d/\n%s\n", c->lon->header,
      count, seq);

  // print to mapping file
  Read* r = c->first;
  fprintf(p->cmap, "%s %s", c->lon->header, r->header);
  for (r = r->next; r != NULL; r = r->next)
    fprintf(p->cmap, "%s%s", COM, r->header);
  fprintf(p->cmap, "\n");

}

/* void makePers()
 * Makes the "consensus" sequence that can be
 *   analyzed by a chimera-checking program.
 */
void makePers(char* seq, char* perSeq, int start) {
  int i;
  for (i = 0; seq[i + start] != '\0'; i++)
    perSeq[i] = seq[i + start];
  perSeq[i] = '\0';
}

/* void makeSeq()
 * Constructs a sequence for an irregular
 *   flow order.
 */
void makeSeq(Read* r, float* flow, int length,
    char* seq, char* order, char* not, int opt) {

  // create sequence
  int start = r->start;
  int len = 0, cons = 0;
  for (int i = 0; i < length; i++) {
    float val = flow[i];

    if (val > MIN) {

      for (; val > MIN; val--)
        seq[len++] = order[i + start];
      cons = 0;
      not[cons++] = order[i + start];

    } else {

      // ambiguous base
      int j;
      for (j = 0; j < cons; j++)
        if (not[j] == order[i + start])
          break;
      if (j == cons) {
        not[cons] = order[i + start];
        if (++cons == NUC && (len || !opt)) {
          seq[len++] = 'N';
          cons = 0;
          not[cons++] = order[i + start];
        }
      }

    }

  }
  seq[len] = '\0';
}

/* void makeIrreg(FILE*)
 * Same as makeReads, but for irregular flow order.
 */
void makeIrreg(FILE* out, Primer* p, Cluster* c, char* seq,
    char* perSeq, char* order, char* trim, int opt,
    int perOpt, int last, int clus, int print) {

  // adjust flowgram, primer
  int begin;
  int pos = strlen(trim);
  if (opt) {
    c->flows[0] = last > c->flows[0] ? 0.0f :
      c->flows[0] - last;
    begin = 0;
  } else {
    float v = c->flows[0];
    int i;
    for (i = 0; i < last && v > MIN; i++)
      v--;
    begin = i;

    if (i != last) {
      int j;
      for (j = 0; j < last - i; j++)
        trim[pos + j] = p->seq[pos + j];
      trim[pos + j] = '\0';
    }
  }

  char* not = (char*) memalloc(NUC); // checking for Ns

  int count = 0;
  Read* r = c->first;
  while (r != NULL) {
    count++;
    makeSeq(r, c->flows, r->length, seq, order, not, opt);

    // print output
    if (opt)
      fprintf(out, ">%s%s%d %s %s\n%s\n", r->mid->name,
        PER, print + count, r->header, p->name, seq);
    else
      fprintf(out, ">%s\n%s%s%s\n", r->header,
        r->mid->seq, trim, seq);

    // save Perseus sequence
    if (r == c->lon && p->per != NULL)
      makePers(seq, perSeq, begin);

    r = r->next;
  }

  free(not);

  if (p->per != NULL)
    printPers(p, c, perSeq, count, clus, perOpt);

  if (!opt)
    trim[pos] = '\0';
}

/* void makeReads(FILE*)
 * Reconstitutes the reads from the consensus flowgrams
 *   and mapping information. Adds on the mid tag and
 *   primer sequences.
 */
void makeReads(FILE* out, Primer* p, Cluster* c, char* seq,
    char* perSeq, char* order, char* trim, int offset,
    int opt, int perOpt, int last, int clus, int print) {

  // adjust flowgram, primer
  int begin;
  int pos = strlen(trim);
  if (opt) {
    c->flows[0] = last > c->flows[0] ? 0.0f :
      c->flows[0] - last;
    begin = 0;
  } else {
    float v = c->flows[0];
    int i;
    for (i = 0; i < last && v > MIN; i++)
      v--;
    begin = i;

    if (i != last) {
      int j;
      for (j = 0; j < last - i; j++)
        trim[pos + j] = p->seq[pos + j];
      trim[pos + j] = '\0';
    }
  }

  int len = 0, cons = 0, count = 0;

  // build seq for each cluster (reads must be sorted)
  Read* r = c->first;
  for (int i = 0; c->flows[i] != END && r != NULL; i++) {

    // ambiguous base
    cons = c->flows[i] <= MIN ? cons + 1 : 0;
    if (cons == 3 && (len || !opt)) {
      seq[len++] = 'N';
      cons = 0;
    }

    for (; c->flows[i] > MIN; c->flows[i]--)
      seq[len++] = order[i % NUC + offset];

    while (r != NULL && r->length == i + 1) {
      seq[len] = '\0';
      count++;

      // print sequence
      if (opt)
        fprintf(out, ">%s%s%d %s %s\n%s\n", r->mid->name,
          PER, print + count, r->header, p->name, seq);
      else
        fprintf(out, ">%s\n%s%s%s\n", r->header,
          r->mid->seq, trim, seq);

      // save Perseus sequence
      if (r == c->lon && p->per != NULL)
        makePers(seq, perSeq, begin);

      r = r->next;
    }
  }

  if (p->per != NULL)
    printPers(p, c, perSeq, count, clus, perOpt);

  if (!opt)
    trim[pos] = '\0';
}

/* int findOffset(char*)
 * Returns an int indicating the last base of the
 *   primer. Trims the last homopolymer from the
 *   primer (returned as input char*).
 */
int findOffset(char* primer, char* trim, char* order) {
  int i = strlen(primer) - 1;
  char last = primer[i];
  int offset = -1;
  for (int j = 0; j < NUC; j++)
    if (last == order[j]) {
      offset = j;
      break;
    }
  if (offset == -1)
    exit(error("", ERRPRIM));

  for (; i && primer[i - 1] == last; i--) ;
  trim[i] = '\0';
  for (--i; i > -1; i--) {
    if (primer[i] == 'A' || primer[i] == 'C' ||
       primer[i] == 'G' || primer[i] == 'T')
      trim[i] = primer[i];

    // ambiguous bases: arbitrary preference for A, then C, G
    else if (primer[i] == 'S' || primer[i] == 'Y' ||
        primer[i] == 'B')
      trim[i] = 'C';
    else if (primer[i] == 'K')
      trim[i] = 'G';
    else
      trim[i] = 'A';
  }
  return offset;
}

/* void printPrimer()
 * Prints the denoised fasta files, consensus
 *   flowgrams and mapping files.
 */
void printPrimer(FILE* out, Primer* p, int len,
    char* order, int denpOpt, int midOpt, int perOpt,
    int reg, int* tclus, int* tread, int statusOpt) {

  char* trim = (char*) memalloc(MIDPRIM);
  int offset = findOffset(p->seq, trim, order);
  int last = strlen(p->seq) - strlen(trim);

  char* seq = (char*) memalloc(len);
  char* perSeq = NULL;
  if (p->per != NULL)
    perSeq = (char*) memalloc(len);
  Read* dummy = (Read*) memalloc(sizeof(Read));
  int clus = 0, read = 0;

  for (Cluster* c = p->head; c != NULL; c = c->next) {

    // skip empty clusters
    if (c->first == NULL)
      continue;

    // sort reads by size
    mergeSort2(c, dummy);

    // print outputs
    if (denpOpt) {

      // print to flow file
      for (int i = 0; c->flows[i] != END; i++)
        fprintf(p->den, "%f ", c->flows[i]);
      fprintf(p->den, "\n");

      // print to mapping file
      for (Read* r = c->first; r != NULL; r = r->next)
        fprintf(p->map, "%s%s%s%s%d ", r->header, SEP,
          r->mid->name, SEP, r->length);
      fprintf(p->map, "\n");

    }

    // print to fasta file
    if (reg)
      makeReads(out, p, c, seq, perSeq, order, trim,
        offset, midOpt, perOpt, last, clus,
        *tread + read);
    else
      makeIrreg(out, p, c, seq, perSeq, order, trim,
        midOpt, perOpt, last, clus, *tread + read);

    read += clusSize2(c);
    clus++;
  }

  if (statusOpt)
    printf("\n%s%sReads: %s%10d\n%s%sClusters:%10d",
      TAB, TAB, TAB, read, TAB, TAB, clus);
  *tclus += clus;
  *tread += read;

  free(dummy);
  free(trim);
  free(seq);
  if (p->per != NULL)
    free(perSeq);
}

/* Read* makeRead(char*, int)
 * Creates a Read in memory.
 */
Read* makeRead(char* header, Midtag* m, int start,
    int last, int numFlows, int reg, int trieOpt) {
  Read* r = (Read*) memalloc(sizeof(Read));
  r->header = (char*) memalloc(1 + strlen(header));
  strcpy(r->header, header);
  r->next = NULL;
  r->mid = m;
  if (!trieOpt) {
    r->length = last - start;
    r->flow = (float*) memalloc(numFlows * sizeof(float));
  }
  if (!reg)
    r->start = start;
  return r;
}

/* void adjustFlow(float*, int, int)
 * Removes the 5' end (key, mid tag, and primer)
 *   from the flowgram.
 */
void adjustFlow(float* flow, float* copy, int start,
    int last) {
  int i;
  for (i = 0; i < last - start; i++)
    copy[i] = flow[i + start];
  copy[i] = END;
}

/* void addOn(Read*, Cluster*, int)
 * Adds flows from the read to the cluster
 *   (for second iteration).
 */
void addOn(Read* r, Cluster* c, int j) {
  for (; r->flow[j] != END; j++)
    c->flows[j] = r->flow[j];
  c->flows[j] = END;
}

/* Cluster* checkClusters(Primer*, float*)
 * Checks the clusters of the given Primer* for
 *   a match of the given flowgram.
 */
Cluster* checkClusters(Primer* p, Read* r,
    float** inter, FILE* missFile, int** miss,
    int iter) {
  Cluster* c = p->head;
  while (c != NULL) {

    for (int j = 0; ; j++) {

      if (r->flow[j] == END || c->flows[j] == END) {
        if (c->flows[j] == END && iter)
          addOn(r, c, j);
        return c;
      }

      int i = (c->flows[j] + 0.005f) * 100.0f;
      if (r->flow[j] < inter[i][0] || r->flow[j] > inter[i][1]) {
        if (iter && missFile != NULL)
          miss[i][(int) ((r->flow[j] + 0.005f) * 100.0f)]++;
        break;
      }
    }

    c = c->next;
  }
  return NULL;
}

/* void makeCluster(Primer*, Read*, float*, int)
 * Creates a new cluster from the given read and flowgram.
 */
void makeCluster(Primer* p, Read* r, int numFlows, int iter) {
  Cluster* neo = (Cluster*) memalloc(sizeof(Cluster));
  if (iter) {
    neo->first = neo->lon = r;
    r->next = NULL;
  } else
    neo->first = neo->lon = NULL;

  // make weights
  neo->weight = (int*) memalloc(numFlows * sizeof(int));
  neo->flows = (float*) memalloc(numFlows * sizeof(float));
  int k;
  for (k = 0; r->flow[k] != END; k++) {
    neo->flows[k] = r->flow[k];
    neo->weight[k] = 1;
  }
  neo->weight[k] = 0;
  neo->flows[k] = END;

  neo->next = NULL;
  if (p->head == NULL)
    p->head = neo;
  else
    p->tail->next = neo;
  p->tail = neo;
}

/* void addToCluster(Read*, Cluster*)
 * Adds the given read to the given cluster by making
 *   a weighted average of the flow values.
 */
void addToCluster(Read* r, Cluster* c, int iter) {
  if (iter) {
    if (c->lon == NULL || r->length > c->lon->length)
      c->lon = r;
    r->next = c->first;
    c->first = r;
  } else {
    int i;
    for (i = 0; i < r->length && c->weight[i]; i++) {
      c->flows[i] = (c->weight[i] * c->flows[i] +
        r->flow[i]) / (c->weight[i] + 1);
      c->weight[i]++;
    }

    if (!c->weight[i]) {
      for (; i < r->length; i++) {
        c->flows[i] = r->flow[i];
        c->weight[i] = 1;
      }
      c->flows[i] = END;
      c->weight[i] = 0;
    }
  }
}

/* void denRead()
 * Denoises a given read.
 */
void denRead(Primer* p, Read* r, float** inter,
    FILE* missFile, int** miss, int numFlows,
    int iter) {
  Cluster* c = checkClusters(p, r, inter,
    missFile, miss, iter);
  if (c == NULL)
    makeCluster(p, r, numFlows, iter);
  else
    addToCluster(r, c, iter);
}

/* Read* denoise(char*, Midtag*, float*)
 * Controls the first iteration of the denoising
 *   process (for the combined clean-denoise analysis).
 */
Read* denoise(char* header, Midtag* m, float* flow,
    int start, int last, int numFlows, float** inter,
    int reg, FILE* missFile, int** miss) {

  Read* r = makeRead(header, m, start, last,
    numFlows, reg, 0);
  adjustFlow(flow, r->flow, start, last);

  denRead(m->prim, r, inter, missFile, miss,
    numFlows, 0);
  return r;
}

/* int checkOrder(char*)
 * Checks if flow order is regular (repetitive).
 */
int checkOrder(char* order, int numFlows) {

  int len = strlen(order);
  if (len < NUC * 2)
    exit(error("", ERRORDER));

  for (int i = 0; i < NUC; i++) {
    if (order[i] != 'A' && order[i] != 'C' &&
        order[i] != 'G' && order[i] != 'T')
      exit(error("", ERRORDER));
    for (int j = NUC; j < len - i; j += NUC)
      if (order[i] != order[i + j]) {
        if (len < numFlows)
          exit(error("", ERRORDER));
        else
          return 0;
      }
  }

  order[NUC * 2] = '\0'; // truncate if repetitive
  return 1;
}

/* Midtag* detMid(Primer*, char*)
 * Parses the Read's header to get the mid tag
 *   and starting flow.
 */
Midtag* detMid(Primer* p, Read* r, int reg) {
  char* name = strtok(r->header, SEP);
  if (name == NULL)
    exit(error("", ERRHEAD));

  // get midtag name
  char* end = reg ? DELIM : SEP;
  name = strtok(NULL, end);
  if (name == NULL)
    exit(error("", ERRHEAD));

  // find midtag
  Midtag* m = p->first;
  while (m != NULL && strcmp(m->name, name))
    m = m->next;
  if (m == NULL)
    exit(error(name, ERRMID));

  // determine start
  if (!reg) {
    char* start = strtok(NULL, DELIM);
    if (start == NULL)
      exit(error("", ERRHEAD));
    r->start = getInt(start);
  }

  return m;
}

/* int getHeader(FILE*)
 * Gets the header information (number of flows,
 *   flow order) from the cleaned flowgram file.
 */
int getHeader(FILE* file, char** order) {
  if (fgets(line, MAX_SIZE, file) == NULL)
    exit(error("", ERRNFLOWS));
  char* in = strtok(line, DELIM);
  if (in == NULL)
    exit(error("", ERRNFLOWS));

  int flows = getInt(in);
  if (!flows)
    exit(error("", ERRNFLOWS));

  *order = (char*) memalloc(flows + 1);
  in = strtok(NULL, DELIM);
  if (in == NULL)
    exit(error("", ERRNFLOWS));
  strcpy(*order, in);

  return flows;
}

/* int getFlows(float*)
 * Loads the flow values from the input line.
 *   Returns the number of good flow values.
 */
int getFlows(float* flow, float max) {
  char* in = strtok(NULL, DELIM);
  int count = 0;
  while (in != NULL) {
    flow[count] = getFloat(in);
    if (flow[count] > max)
      flow[count] = max;
    count++;
    in = strtok(NULL, DELIM);
  }
  flow[count] = END;

  return count;
}

/* Read* getRead()
 * Initializes the next Read from already
 *   cleaned flowgrams.
 */
Read* getRead(Primer* p, int numFlows, int reg,
    float* flow, float max) {

  if (fgets(line, MAX_SIZE, p->out) == NULL)
    return NULL;
  char* head = strtok(line, DELIM);
  if (head == NULL)
    exit(error("", ERRHEAD));

  Read* r = (Read*) memalloc(sizeof(Read));
  r->next = NULL;
  r->header = (char*) memalloc(1 + strlen(head));
  strcpy(r->header, head);

  if (flow == NULL) {
    r->flow = (float*) memalloc(numFlows * sizeof(float));
    r->length = getFlows(r->flow, max);
  } else
    r->length = getFlows(flow, max);

  r->mid = detMid(p, r, reg);

  return r;
}

/* void secondIter()
 * Controls the second iteration of denoising.
 */
void secondIter(Primer* p, float** inter,
    FILE* missFile, int** miss, int numFlows,
    int total) {

  if (total)
    printf("\n%sPrimer %s (second iteration)\n",
      TAB, p->name);

  // sort clusters by size
  Cluster* dummy = (Cluster*) memalloc(sizeof(Cluster));
  mergeSort(p, dummy);
  free(dummy);

  // denoise reads -- 2nd iteration
  int count = 0;
  Read* r = p->dummy->next;
  Read* next;
  while (r != NULL) {
    if (total)
      printf("%s%sCompleted %.1f%%\r", TAB, TAB,
        ++count * 100.0f / total);
    next = r->next;
    denRead(p, r, inter, missFile, miss, numFlows, 1);
    free(r->flow);
    r = next;
  }

}

/* void printDenoise()
 * Controls the second iteration and output
 *   for the combined filter-denoise analysis.
 */
void printDenoise(FILE* out, char* order, int denpOpt,
    int midOpt, int perOpt, int reg, FILE* missFile,
    int** miss, float** inter, int numFlows,
    float max, int statusOpt) {

  if (statusOpt)
    printf("Denoising:");
  int tclus = 0, tread = 0;

  for (Primer* p = primo; p != NULL; p = p->next) {
    int total = 0;
    if (statusOpt)
      for (Read* r = p->dummy->next; r != NULL; r = r->next)
        total++;
    secondIter(p, inter, missFile, miss, numFlows, total);
    printPrimer(out, p, 2 * numFlows, order, denpOpt,
      midOpt, perOpt, reg, &tclus, &tread, statusOpt);

    freeClus(p->head);
    p->head = NULL;
  }

  if (statusOpt)
    printf("\n%sTotal\n%s%sReads: %s%10d\n%s%sClusters:%10d\n",
      TAB, TAB, TAB, TAB, tread, TAB, TAB, tclus);

  // print misses array
  if (missFile != NULL) {
    int maxm = (int) ((max + 0.01f) * 100.0f);
    printMiss(missFile, miss, maxm);
  }

  if (fclose(out))
    exit(error("", ERRCLOSE));
}

/* void printTrPers()
 * Prints the output that can be analyzed for
 *   chimeras.
 */
void printTrPers(Read* r, char* seq, int begin,
    Read** map, int mapPos, int perOpt, int* leaf) {

  // print to fasta file
  if (perOpt)
    fprintf(r->mid->prim->per, ">%s%s%d%s%d\n%s\n",
      r->header, PER, *leaf, PER, mapPos, seq + begin);
  else
    fprintf(r->mid->prim->per, ">%s /ab=%d/\n%s\n",
      r->header, mapPos, seq + begin);

  // print to mapping file
  FILE* out = r->mid->prim->cmap;
  fprintf(out, "%s %s", r->header, map[0]->header);
  for (int i = 1; i < mapPos; i++)
    fprintf(out, "%s%s", COM, map[i]->header);
  fprintf(out, "\n");

}

/* void printFlow()
 * Prints the denoised flowgram for a read.
 */
void printFlow(Read* r, float* flow, int len,
    int reg) {
  FILE* out = r->mid->prim->den;
  fprintf(out, "%s%s%s", r->header, SEP,
    r->mid->name);
  if (!reg)
    fprintf(out, "%s%d", SEP, r->start);

  for (int i = 0; i < len; i++)
    fprintf(out, " %f", flow[i]);
  fprintf(out, "\n");
}

/* void printIrreg()
 * Prints the output for this node
 *   (irregular flow order).
 */
void printIrreg(FILE* out, Node* n, float* flow,
    int flowPos, char* seq, char* order,
    char* trim, int midOpt, int total, int* count,
    int denpOpt, char* not, int begin,
    Read** map, int* mapPos, int perOpt, int* leaf) {

  Read* prev = NULL;
  for (Read* r = n->first; r != NULL; r = r->next) {
    (*count)++;
    prev = r;

    // make sequence
    makeSeq(r, flow, flowPos, seq, order, not, midOpt);

    // print sequence
    if (midOpt)
      fprintf(out, ">%s%s%d %s %s\n%s\n", r->mid->name, PER,
        total + *count, r->header, r->mid->prim->name, seq);
    else
      fprintf(out, ">%s\n%s%s%s\n", r->header,
        r->mid->seq, trim, seq);

    // print flowgram
    if (denpOpt) {
      float save;
      if (midOpt) {
        save = flow[0];
        flow[0] += strlen(r->mid->prim->seq) - strlen(trim);
      }
      printFlow(r, flow, flowPos, 0);
      if (midOpt)
        flow[0] = save;
    }

    if (map != NULL)
      map[(*mapPos)++] = r;
  }

  // print chimera output (if no child)
  if (prev != NULL && n->child == NULL) {
    if (map != NULL)
      printTrPers(prev, seq, begin, map, *mapPos,
        perOpt, leaf);
    (*leaf)++;
  }

}

/* void addFlow()
 * Controls the output for an irregular
 *   flow order.
 */
void addFlow(FILE* out, Node* n, float* flow,
    int flowPos, char* seq, char* order,
    char* trim, int midOpt, int total, int* count,
    int denpOpt, char* not, int begin,
    Read** map, int mapPos, int perOpt, int* leaf) {

  if (n == NULL)
    return;

  int oldFlPos = flowPos;
  int oldMapPos = mapPos;

  for (int i = n->st; i < n->end; i++)
    flow[flowPos++] = n->flow[i];

  printIrreg(out, n, flow, flowPos, seq, order,
    trim, midOpt, total, count, denpOpt, not,
    begin, map, &mapPos, perOpt, leaf);

  addFlow(out, n->child, flow, flowPos, seq, order,
    trim, midOpt, total, count, denpOpt, not,
    begin, map, mapPos, perOpt, leaf);

  addFlow(out, n->next, flow, oldFlPos, seq, order,
    trim, midOpt, total, count, denpOpt, not,
    begin, map, oldMapPos, perOpt, leaf);
}

/* void printNode()
 * Prints the output for this node.
 */
void printNode(FILE* out, Node* n, char* seq, int pos,
    int level, char* trim, int midOpt, int total,
    int* count, float* flow, int begin, Read** map,
    int* mapPos, int perOpt, int* leaf) {

  Read* prev = NULL;
  for (Read* r = n->first; r != NULL; r = r->next) {
    seq[pos] = '\0';
    (*count)++;
    prev = r;

    // print sequence
    if (midOpt)
      fprintf(out, ">%s%s%d %s %s\n%s\n", r->mid->name, PER,
        total + *count, r->header, r->mid->prim->name, seq);
    else
      fprintf(out, ">%s\n%s%s%s\n", r->header,
        r->mid->seq, trim, seq);

    // print flowgram
    if (flow != NULL)
      printFlow(r, flow, level, 1);

    if (map != NULL)
      map[(*mapPos)++] = r;
  }

  // print chimera output (if no child)
  if (prev != NULL && n->child == NULL) {
    if (map != NULL)
      printTrPers(prev, seq, begin, map, *mapPos,
        perOpt, leaf);
    (*leaf)++;
  }

}

/* void buildSeq()
 * Adds to the sequence from the given Node.
 */
void buildSeq(Node* n, char* seq, int* pos, int* level,
    char* order, int offset, int* no, float* flow,
    int opt) {

  for (int i = n->st; i < n->end; i++) {

    float val = n->flow[i];
    if (val > MIN) {
      for (; val > MIN; val--)
        seq[(*pos)++] = order[*level % NUC + offset];
      *no = 0;
    } else if (++*no == 3 && (*pos || !opt)) {
      seq[(*pos)++] = 'N';
      *no = 0;
    }

    if (flow != NULL && *level)
      flow[*level] = n->flow[i];

    (*level)++;
  }

}

/* void addBase()
 * Controls the output at a node.
 */
void addBase(FILE* out, Node* n, char* seq, int pos,
    int level, char* order, int offset, char* trim,
    int midOpt, int total, int* count, int no,
    float* flow, int begin, Read** map, int mapPos,
    int perOpt, int* leaf) {

  if (n == NULL)
    return;

  int oldPos = pos;
  int oldLevel = level;
  int oldNo = no;
  int oldMapPos = mapPos;

  buildSeq(n, seq, &pos, &level, order, offset,
    &no, flow, midOpt);

  printNode(out, n, seq, pos, level, trim, midOpt,
    total, count, flow, begin, map, &mapPos, perOpt,
    leaf);

  addBase(out, n->child, seq, pos, level, order,
    offset, trim, midOpt, total, count, no, flow,
    begin, map, mapPos, perOpt, leaf);

  addBase(out, n->next, seq, oldPos, oldLevel, order,
    offset, trim, midOpt, total, count, oldNo, flow,
    begin, map, oldMapPos, perOpt, leaf);
}

/* int nodeSize()
 * Returns the maximum sub-trie size.
 */
int nodeSize(Primer* p) {
  int ans = 0;
  for (Node* n = p->root->child; n != NULL; n = n->next)
    if (ans < n->num)
      ans = n->num;
  return ans;
}

/* void printTrie()
 * Prints the denoised fasta files, consensus
 *   flowgrams and mapping files.
 */
void printTrie(FILE* out, Primer* p, int len,
    char* order, int denpOpt, int midOpt, int perOpt,
    int reg, int* tclus, int* tread, int statusOpt) {

  char* trim = (char*) memalloc(MIDPRIM);
  int offset = findOffset(p->seq, trim, order);
  int last = strlen(p->seq) - strlen(trim);

  char* seq = (char*) memalloc(len);

  Read** map = NULL;
  if (p->per != NULL)
    map = (Read**) memalloc(nodeSize(p) * sizeof(Read*));

  float* flow = NULL;
  if (!reg || denpOpt)
    flow = (float*) memalloc(len / 2 * sizeof(float));

  char* not;
  if (!reg)
    not = (char*) memalloc(NUC);

  int begin;
  int st = strlen(trim);
  int count = 0, leaf = 0;

  Node* n = p->root->child;
  while (n != NULL) {

    if (denpOpt)
      flow[0] = n->flow[n->st];

    // adjust flowgram, primer
    if (midOpt) {
      n->flow[n->st] = last > n->flow[n->st] ? 0.0f :
        n->flow[n->st] - last;
      begin = 0;
    } else {
      float v = n->flow[n->st];
      int i;
      for (i = 0; i < last && v > MIN; i++)
        v--;
      begin = i;

      if (i != last) {
        int j;
        for (j = 0; j < last - i; j++)
          trim[st + j] = p->seq[st + j];
        trim[st + j] = '\0';
      }
    }

    // analyze trie
    int pos = 0, level = 0, no = 0, mapPos = 0;

    if (!reg) {

      // start flowgram
      for (int i = n->st; i < n->end; i++)
        flow[pos++] = n->flow[i];

      // print output (if necessary)
      printIrreg(out, n, flow, pos, seq, order,
        trim, midOpt, *tread, &count, denpOpt,
        not, begin, map, &mapPos, perOpt, &leaf);

      // enter trie
      addFlow(out, n->child, flow, pos, seq,
        order, trim, midOpt, *tread, &count,
        denpOpt, not, begin, map, mapPos,
        perOpt, &leaf);

    } else {

      // start sequence
      buildSeq(n, seq, &pos, &level, order, offset,
        &no, flow, midOpt);

      // print output (if necessary)
      printNode(out, n, seq, pos, level, trim, midOpt,
        *tread, &count, flow, begin, map, &mapPos, perOpt,
        &leaf);

      // enter trie
      addBase(out, n->child, seq, pos, level, order,
        offset, trim, midOpt, *tread, &count, no,
        flow, begin, map, mapPos, perOpt, &leaf);

    }

    if (!midOpt)
      trim[st] = '\0';

    n = n->next;
  }

  if (statusOpt)
    printf("\n%s%sReads: %s%s%10d\n%s%sLeaf nodes:%10d",
      TAB, TAB, TAB, TAB, count, TAB, TAB, leaf);
  *tclus += leaf;
  *tread += count;

  free(trim);
  free(seq);
  if (p->per != NULL)
    free(map);
  if (flow != NULL)
    free(flow);
  if (!reg)
    free(not);
}

/* void printNodes()
 * Dumps the trie (for debugging).
 */
void printNodes(Node* first, int level) {
  for (Node* n = first; n != NULL; n = n->next) {
    for (int i = 0; i < level; i++)
      printf("%s", TAB);
    printf("%.3f to %.3f (%d nodes)\n", n->flow[n->st],
      n->flow[n->end - 1], n->end - n->st);
    printNodes(n->child, level + 1);
  }
}

/* float* nodeFlow()
 * Copies a subset of flows for a node.
 */
float* nodeFlow(float* flow, int pos, int size) {
  float* ans = (float*) memalloc(size * sizeof(float));
  for (int i = 0; i < size; i++)
    ans[i] = flow[i + pos];
  return ans;
}

/* void makeNode()
 * Makes a node.
 */
Node* makeNode(Node* par, float* flow, int pos,
    int size) {
  Node* n = (Node*) memalloc(sizeof(Node));
  n->num = 1;
  n->first = NULL;
  n->next = par->child;
  par->child = n;
  n->child = NULL;

  n->flow = nodeFlow(flow, pos, size);
  n->end = size;
  n->st = 0;

  return n;
}

/* Node* matchFlow()
 * Finds the best match of a list of (children)
 *   nodes to the query.
 */
Node* matchFlow(Node* n, float query, float** inter,
    FILE* missFile, int** miss) {

  Node* ans = NULL;
  float dist = DEFMAXFLOW;

  while (n != NULL) {
    float val = n->flow[n->st];
    int i = (val + 0.005f) * 100.0f;

    if (query < inter[i][0] || query > inter[i][1] ||
        query - val >= dist || query - val <= -dist) {

      if (missFile != NULL)
        miss[i][(int) ((query + 0.005f) * 100.0f)]++;

    } else {
      ans = n;
      float d = query - val;
      dist = d < 0 ? -d : d;
    }

    n = n->next;
  }

  return ans;
}

/* void splitNode()
 * Splits a node by creating a new child.
 */
void splitNode(Node* par, int size) {

  Node* n = (Node*) memalloc(sizeof(Node));
  n->num = par->num;
  n->first = par->first;
  n->next = NULL;
  n->child = par->child;

  // alloc new float*
  if (size < (par->end - par->st) / 2.0) {
    n->flow = par->flow;
    n->st = par->st + size;
    n->end = par->end;

    // parent gets new float*
    par->flow = nodeFlow(par->flow, par->st, size);
    par->st = 0;
    par->end = size;

  } else {

    // child gets new float*
    n->st = 0;
    n->end = par->end - par->st - size;
    n->flow = nodeFlow(par->flow, par->st + size, n->end);
    par->end = par->st + size;

  }

  par->first = NULL;
  par->child = n;
  par->num++;

}

// declaration for recursive call from moreFlow()
Node* checkNode(Node*, float*, int, int, float**,
  FILE*, int**);

/* Node* moreFlow()
 * Checks remaining flows of node.
 */
Node* moreFlow(Node* n, float* flow, int pos, int len,
    float** inter, FILE* missFile, int** miss) {

  int j;
  for (j = 1; j < n->end - n->st; j++) {

    // if query is short, split node
    if (pos + j == len) {
      splitNode(n, j);
      return n;
    }

    int i = (n->flow[n->st + j] + 0.005f) * 100.0f;

    if (flow[pos + j] < inter[i][0] ||
        flow[pos + j] > inter[i][1]) {

      if (missFile != NULL)
        miss[i][(int) ((flow[pos + j] + 0.005f) * 100.0f)]++;

      // split node, and make new child
      splitNode(n, j);
      return makeNode(n, flow, pos + j, len - pos - j);

    } else

      // make weighted average, and continue
      n->flow[n->st + j] = (n->flow[n->st + j] * n->num +
        flow[pos + j]) / (n->num + 1);

  }

  // still need to add 1 to n->num
  n->num++;

  // recurse on child node
  return checkNode(n, flow, pos + j, len, inter,
    missFile, miss);
}

/* Node* checkNode()
 * Places a read into the trie.
 */
Node* checkNode(Node* n, float* flow, int pos, int len,
    float** inter, FILE* missFile, int** miss) {

  // base case: end of flowgram
  if (pos == len)
    return n;

  // check for matching node
  Node* m = matchFlow(n->child, flow[pos], inter,
    missFile, miss);

  // if no match, make a new node
  if (m == NULL)
    return makeNode(n, flow, pos, len - pos);

  // make weighted average, check remaining flows
  m->flow[m->st] = (m->flow[m->st] * m->num + flow[pos]) /
    (m->num + 1);

  return moreFlow(m, flow, pos, len, inter,
    missFile, miss);
}

/* void denTrie()
 * Denoises the given read using a trie.
 */
void denTrie(Primer* p, Read* r, int len, float* flow,
    float** inter, FILE* missFile, int** miss) {

  // enter trie
  Node* n = checkNode(p->root, flow, 0, len, inter,
    missFile, miss);

  // add read to node
  r->next = n->first;
  n->first = r;

}

/* void readClean()
 * Controls the denoising process (for the
 *   denoising only option).
 */
void readClean(int midOpt, float** inter,
    int denpOpt, FILE* denFasta, FILE* missFile,
    int** miss, int trieOpt, int perOpt, float max,
    int statusOpt) {

  if (statusOpt)
    printf("Denoising%s:", trieOpt ?
      " (by trie)" : "");
  int tclus = 0, tread = 0;
  float* flow = NULL;

  Primer* p = primo;
  while (p != NULL) {

    int count = 0, total = 0;
    if (statusOpt) {
      printf("\n%sPrimer %s %s\n", TAB, p->name,
        trieOpt ? "" : "(first iteration)");
      while (fgets(line, MAX_SIZE, p->out) != NULL)
        total++;
      total--;
      rewind(p->out);
    }

    // get header information
    char* order = NULL; // flow order
    int numFlows = getHeader(p->out, &order);
    int reg = checkOrder(order, numFlows);
    numFlows++;

    if (trieOpt) {

      flow = (float*) memalloc(numFlows * sizeof(float));

      Read* r = getRead(p, numFlows, reg, flow, max);
      while (r != NULL) {
        if (total)
          printf("%s%sCompleted %.1f%%\r", TAB, TAB,
            ++count * 100.0f / total);
        denTrie(p, r, r->length, flow, inter,
          missFile, miss);
        r = getRead(p, numFlows, reg, flow, max);
      }

      printTrie(denFasta, p, 2 * numFlows, order,
        denpOpt, midOpt, perOpt, reg, &tclus, &tread,
        statusOpt);

      freeTrie(p->root);
      p->root = NULL;
      free(flow);

    } else {

      // denoise reads -- 1st iteration
      Read* r = getRead(p, numFlows, reg, flow, max);
      p->dummy->next = r;
      p->prev = r;
      while (r != NULL) {
        if (total)
          printf("%s%sCompleted %.1f%%\r", TAB, TAB,
            ++count * 100.0f / total);
        denRead(p, r, inter, missFile, miss,
          numFlows, 0);
        r = getRead(p, numFlows, reg, flow, max);
        p->prev->next = r;
        p->prev = r;
      }

      // 2nd iteration
      secondIter(p, inter, missFile, miss, numFlows, total);

      // produce outputs (flowgram, mapping, fasta)
      printPrimer(denFasta, p, 2 * numFlows, order, denpOpt,
        midOpt, perOpt, reg, &tclus, &tread, statusOpt);

      freeClus(p->head);
      p->head = NULL;
    }

    free(order);
    p = p->next;
  }

  if (statusOpt)
    printf("\n%sTotal\n%s%sReads: %s%s%10d\n%s%s%s:%10d\n",
      TAB, TAB, TAB, TAB, trieOpt ? TAB : "", tread, TAB, TAB,
      trieOpt ? "Leaf nodes" : "Clusters", tclus);

  // free inter array
  int maxm = (int) ((max + 0.01f) * 100.0f);
  for (int i = 0; i < maxm; i++)
    free(inter[i]);
  free(inter);

  // print misses array
  if (missFile != NULL)
    printMiss(missFile, miss, maxm);

  if (fclose(denFasta))
    exit(error("", ERRCLOSE));
}

/* void printDenTrie()
 * Prints the output fasta from the combined
 *   filter-denoise analysis.
 */
void printDenTrie(FILE* out, int len, char* order,
    int denpOpt, int midOpt, int perOpt, int reg,
    FILE* missFile, int** miss, float max,
    int statusOpt) {

  if (statusOpt)
    printf("Denoising results:");

  int tclus = 0, tread = 0;
  for (Primer* p = primo; p != NULL; p = p->next) {
    if (statusOpt)
      printf("\n%sPrimer %s", TAB, p->name);
    printTrie(out, p, len, order, denpOpt, midOpt,
      perOpt, reg, &tclus, &tread, statusOpt);

    freeTrie(p->root);
    p->root = NULL;
  }

  if (statusOpt)
    printf("\n%sTotal\n%s%sReads: %s%s%10d\n%s%sLeaf nodes:%10d\n",
      TAB, TAB, TAB, TAB, TAB, tread, TAB, TAB, tclus);

  // print misses array
  if (missFile != NULL) {
    int maxm = (int) ((max + 0.01f) * 100.0f);
    printMiss(missFile, miss, maxm);
  }

  if (fclose(out))
    exit(error("", ERRCLOSE));
}

/* void trieStart()
 * Denoises a read using a trie for the
 *   combined filter-denoise analysis.
 */
void trieStart(char* header, Midtag* m, float* flow,
    int start, int last, float** inter, int reg,
    FILE* missFile, int** miss) {

  Read* r = makeRead(header, m, start, last, 0, reg, 1);
  *(flow + last) = END; // add tag to end

  denTrie(m->prim, r, last - start, flow + start,
    inter, missFile, miss);

}

/* char* getData(FILE*, char*)
 * Returns the value indicated by the given char*.
 */
char* getData(FILE* file, char* lin, char* data) {
  while (fgets(lin, MAX_SIZE, file) != NULL) {
    if (lin[0] == '>')
      break;
    char* label = strtok(lin, COL);
    char* value = strtok(NULL, DELIM);
    if (!strcmp(label, data))
      return value;
  }
  exit(error(data, ERRPARAM));
}

/* int ambig(char, char)
 * Checks ambiguous DNA bases.
 */
int ambig(char x, char y) {
  if (x == 'N' ||
      (x == 'W' && (y == 'A' || y == 'T')) ||
      (x == 'S' && (y == 'C' || y == 'G')) ||
      (x == 'M' && (y == 'A' || y == 'C')) ||
      (x == 'K' && (y == 'G' || y == 'T')) ||
      (x == 'R' && (y == 'A' || y == 'G')) ||
      (x == 'Y' && (y == 'C' || y == 'T')) ||
      (x == 'B' && (y == 'C' || y == 'G' || y == 'T')) ||
      (x == 'D' && (y == 'A' || y == 'G' || y == 'T')) ||
      (x == 'H' && (y == 'A' || y == 'C' || y == 'T')) ||
      (x == 'V' && (y == 'A' || y == 'C' || y == 'G')))
    return 0;
  return 1;
}

/* Midtag* findMid(char*)
 * Finds a mid tag - primer match to the given sequence.
 */
Midtag* findMid(char* seq, int midMis, int primMis) {
  for (Primer* p = primo; p != NULL; p = p->next) {
    for (Midtag* m = p->first; m != NULL; m = m->next) {

      // check mid tag
      int misAllow = midMis;
      int i;
      for (i = 0; m->seq[i] != '\0'; i++)
        if (seq[i] == '\0' || (m->seq[i] != seq[i] &&
            (m->seq[i] == 'A' || m->seq[i] == 'C' ||
            m->seq[i] == 'G' || m->seq[i] == 'T' ||
            ambig(m->seq[i], seq[i])) && --misAllow < 0))
          break;

      // check primer
      if (m->seq[i] == '\0') {
        misAllow = primMis;
        int j;
        for (j = 0; p->seq[j] != '\0'; j++)
          if (seq[i + j] == '\0' || (p->seq[j] != seq[i + j] &&
              (p->seq[j] == 'A' || p->seq[j] == 'C' ||
              p->seq[j] == 'G' || p->seq[j] == 'T' ||
              ambig(p->seq[j], seq[i + j])) &&
              (p->seq[j + 1] == '\0' || --misAllow < 0)))
            break;

        if (p->seq[j] == '\0')
          return seq[i + j] != '\0' ? m : NULL;
      }
    }
  }
  return NULL;
}

/* void loadFlow(char*, float*)
 * Loads the flow values from the given char*
 *   to the given float*.
 */
void loadFlow(char* line, float* flow, int last,
    float max) {
  char* in = strtok(line, CSV);
  int i = 0;
  while (in != NULL && i < last) {
    flow[i] = getFloat(in);
    if (flow[i] > max)
      flow[i] = max;
    i++;
    in = strtok(NULL, CSV);
  }
}

/* int loadInd(char*, int*, int)
 * Loads the flow indexes from the given char*
 *   to the given int*. (Also used to load the
 *   quality scores)
 */
int loadInd(char* line, int* index, int clipl, int len, int end) {
  char* in = strtok(line, CSV);
  int ans;
  for (int i = 0; i < clipl; i++) {
    if (in == NULL)
      exit(error("", ERRLOAD));
    ans = getInt(in);
    in = strtok(NULL, CSV);
  }
  if (getInt(in) == ans)
    ans--;

  int i = 0;
  while (in != NULL && i < len) {
    index[i++] = getInt(in);
    in = strtok(NULL, CSV);
  }
  index[i] = end;  // tag on end
  return ans;
}

/* int checkRev()
 * Checks the seq for a match of the reverse primer.
 */
int checkRev(char* rev, char* seq, int i, int mis) {
  for (int j = 1; rev[j] != '\0'; j++)
    if (seq[i + j] == '\0' || (rev[j] != seq[i + j] &&
        (rev[j] == 'A' || rev[j] == 'C' ||
        rev[j] == 'G' || rev[j] == 'T' ||
        ambig(rev[j], seq[i + j])) && --mis < 0))
      return 0;
  return 1;
}

/* int checkSeq()
 * Checks the read for length, ambiguous bases,
 *   homopolymer runs, reverse primer.
 */
int checkSeq(char* seq, int length, int minLength,
    int maxLength, int maxTrunc, int maxAmbig,
    int okAmbig, int maxHomo, int okHomo, int revmOpt,
    int revqOpt, char* rev, int revMis, int* errNum) {

  // check length
  if (length < minLength)
    *errNum = EMINSLEN;
  else if (length > maxLength)
    *errNum = EMAXSLEN;
  if (*errNum)
    return 0;

  int i, ans = 0, cons = 1;
  for (i = 0; i < length; i++) {

    // truncation length
    if (i == maxTrunc) {
      if (!ans)
        *errNum = EMAXTRLEN;
      break;
    }

    // ambiguous base
    if (seq[i] == 'N') {
      if (--maxAmbig < 0) {
        *errNum = EMAXAMBIG;
        return 0;
      }
      if (--okAmbig < 0 && (i < minLength || !ans)) {
        *errNum = EOKAMBIG;
        if (i < minLength)
          return 0;
        ans = i;
      }
    }

    // homopolymer length
    if (maxHomo || okHomo) {
      cons = (i && seq[i] == seq[i - 1]) ? cons + 1 : 1;
      if (maxHomo && cons == maxHomo + 1) {
        *errNum = EMAXHOMO;
        return 0;
      }
      if (okHomo && cons == okHomo + 1 &&
          (i - okHomo < minLength || !ans)) {
        *errNum = EOKHOMO;
        if (i - okHomo < minLength)
          return 0;
        ans = i - okHomo;
      }
    }

    // reverse primer
    if ((revmOpt || revqOpt) && !ans && seq[i] == rev[0]
        && checkRev(rev, seq, i, revMis)) {
      *errNum = EREVERSE;
      if (i < minLength)
        return 0;
      ans = i;
    }

  }

  if (revqOpt && *errNum != EREVERSE) {
    *errNum = EREVERSE;
    return 0;
  }

  return ans ? ans : i;
}

/* int checkQual()
 * Checks quality scores for a read.
 */
int checkQual(int* qual, int len, int minLength,
    float averageQual, int window, float windowAvg,
    int opt, int* errNum) {

  int i, ans = 0, total = 0;
  for (i = 0; i < len; i++) {
    if (window && i + window - 1 < len) {
      int wind = 0;
      for (int j = 0; j < window; j++)
        wind += qual[i + j];
      if (!ans && (float) wind / window < windowAvg) {
        *errNum = EWINDOW;
        ans = i < minLength || opt ? -1 : i;
      }
    }
    total += qual[i];
  }

  if ((float) total / len < averageQual) {
    *errNum = EAVGQUAL;
    return 0;
  }
  if (ans == -1)
    return 0;

  return ans ? ans : i;
}

/* int checkFlows()
 * Checks flows for a read.
 */
int checkFlows(float* flowgram, int start, int end,
    int minFlows, int maxFlows, float minNoise, float maxNoise,
    float maxValue, float noFlow, int* errNumF) {

  if (end < minFlows) {
    *errNumF = EMINFLEN;
    return 0;
  }

  int i, cons = 0;
  for (i = start; i < end; i++) {
    if (i == maxFlows) {
      *errNumF = EMAXFLEN;
      return i;
    }
    if (maxNoise && flowgram[i] >= minNoise && flowgram[i] <= maxNoise) {
      *errNumF = EFLOWINT;
      return i < minFlows ? 0 : i;
    }
    if (flowgram[i] > maxValue) {
      *errNumF = EMAXNVAL;
      return i < minFlows ? 0 : i;
    }

    if (noFlow) {
      cons = flowgram[i] < noFlow ? cons + 1 : 0;
      if (cons == NUC) {
        *errNumF = ENOFLOW;
        return i - 3 < minFlows ? 0 : i - 3;
      }
    }
  }

  return i;
}

/* void printClean()
 * Prints output from cleaning step.
 */
void printClean(FILE* fasta, char* header, Midtag* mid,
    float* flowgram, int start, int last, char* seq,
    int midOpt, int reg, int count) {

  // print flowgram
  fprintf(mid->prim->out, "%s%s%s", header, SEP, mid->name);
  if (!reg)
    fprintf(mid->prim->out, "%s%d", SEP, start);
  for (int i = start; i < last; i++)
    fprintf(mid->prim->out, " %.2f", flowgram[i]);
  fprintf(mid->prim->out, "\n");

  // print fasta sequence
  if (fasta != NULL) {
    if (midOpt) {
      int add = strlen(mid->seq) + strlen(mid->prim->seq);
      fprintf(fasta, ">%s%s%d %s %s\n%s\n", mid->name, PER,
        count, header, mid->prim->name, seq + add);
    } else
      fprintf(fasta, ">%s\n%s\n", header, seq);
  }
}

/* int naCheck(int, int, int)
 * Checks if category is "n/a".
 */
int naCheck(int i, int j, int opt, int val, int revOpt) {
  if ((i == ELIM && !val && (j == EMAXTRLEN ||
      j == EMAXFLEN || (j == EREVERSE && revOpt))) ||
      (i == TRUNC && (j == EMINSLEN || j == EMAXSLEN ||
      j == EMAXAMBIG || j == EMAXHOMO || j == EAVGQUAL ||
      j == EMINFLEN || (j == EWINDOW && opt))))
    return 1;
  return 0;
}

/* void printSamples()
 * Prints detailed filtering counts for each sample.
 */
void printSamples(FILE* err, int elim[][ETCAT][ERRCAT],
    int match[][ETCAT], char** categ, int opt, int revOpt) {
  fprintf(err, "Breakdown by sample:\n");
  for (Primer* p = primo; p != NULL; p = p->next) {
    fprintf(err, "\n%s\n", p->name);
    for (int i = 0; i < ERRCAT; i++)
      if (categ[i] != NULL)
        fprintf(err, "\t%s", categ[i]);
    fprintf(err, "\t%s", TOTAL);
    int tot[ETCAT][ERRCAT + 2] = {{0}};

    // print totals for each sample
    for (Midtag* m = p->first; m != NULL; m = m->next) {
      fprintf(err, "%s\n", m->name);
      for (int j = ELIM; j < ETCAT; j++) {
        j ? fprintf(err, STRUNC) : fprintf(err, SELIM);
        int k, total = 0;
        for (k = 0; k < ERRCAT; k++)
          if (categ[k] != NULL) {
            if (naCheck(j, k, opt, elim[m->num][j][k], revOpt))
              fprintf(err, "\t%s", NA);
            else {
              fprintf(err, "\t%d", elim[m->num][j][k]);
              total += elim[m->num][j][k];
              tot[j][k] += elim[m->num][j][k];
            }
          }
        fprintf(err, "\t%d", total);
        tot[j][k++] += total;
        if (j) {
          fprintf(err, "\t%s\t%d\n", SPRINT, match[m->num][PRINT]);
          tot[j][k] += match[m->num][PRINT];
        } else {
          fprintf(err, "\t%s\t%d\n", SMATCH, match[m->num][MATCH]);
          tot[j][k] += match[m->num][MATCH];
        }
      }
    }

    // print totals for the primer
    fprintf(err, "\n%s Totals\n", p->name);
    for (int j = ELIM; j < ETCAT; j++) {
      j ? fprintf(err, STRUNC) : fprintf(err, SELIM);
      int k;
      for (k = 0; k < ERRCAT; k++)
        if (categ[k] != NULL)
          naCheck(j, k, opt, tot[j][k], revOpt) ?
            fprintf(err, "\t%s", NA) :
            fprintf(err, "\t%d", tot[j][k]);
      fprintf(err, "\t%d", tot[j][k++]);
      j ? fprintf(err, "\t%s\t%d\n", SPRINT, tot[j][k]) :
        fprintf(err, "\t%s\t%d\n", SMATCH, tot[j][k]);
    }
  }

}

/* void printMatch()
 * Prints the match/print counts.
 */
void printMatch(FILE* err, int count, int match[][ETCAT],
    int statusOpt) {

  if (statusOpt) {
    printf("\n%s%s%s%s%10d\n", TAB, COUNT,
      TAB, TAB, count);
    for (int i = MATCH; i < ETCAT; i++) {
      if (i)
        printf("%s%s %s%s", TAB, SPRINT, TAB, TAB);
      else
        printf("%s%s", TAB, SMATCH);
      int k = 0, total = 0;
      for (Primer* p = primo; p != NULL; p = p->next)
        for (Midtag* m = p->first; m != NULL; m = m->next)
          total += match[k++][i];
      printf("%10d\n", total);
    }
  }

  if (err != NULL) {
    fprintf(err, "\n%s\t%d\n", COUNT, count);
    for (Primer* p = primo; p != NULL; p = p->next)
      fprintf(err, "\t%s", p->name);
    fprintf(err, "\t%s", TOTAL);
    for (int i = MATCH; i < ETCAT; i++) {
      fprintf(err, "%s", i ? SPRINT : SMATCH);
      int k = 0, total = 0;
      for (Primer* p = primo; p != NULL; p = p->next) {
        int tot = 0;
        for (Midtag* m = p->first; m != NULL; m = m->next)
          tot += match[k++][i];
        fprintf(err, "\t%d", tot);
        total += tot;
      }
      fprintf(err, "\t%d\n", total);
    }
    fprintf(err, "\n");
  }
}

/* void printStats()
 * Prints the cleaning data.
 */
void printStats(FILE* err, int elim[][ETCAT][ERRCAT],
    int match[][ETCAT], char** categ, int opt, int count,
    int samples, int revmOpt, int statusOpt) {

  // print elim/trunc counts
  if (err != NULL) {
    for (int i = 0; i < ERRCAT; i++)
      if (categ[i] != NULL)
        fprintf(err, "\t%s", categ[i]);
    fprintf(err, "\t%s", TOTAL);

    for (int i = ELIM; i < ETCAT; i++) {
      i ? fprintf(err, STRUNC) : fprintf(err, SELIM);
      int total = 0;
      for (int j = 0; j < ERRCAT; j++) {
        if (categ[j] != NULL) {
          int subt = 0;
          for (int k = 0; k < samples; k++)
            subt += elim[k][i][j];
          naCheck(i, j, opt, subt, revmOpt) ?
            fprintf(err, "\t%s", NA) :
            fprintf(err, "\t%d", subt);
          total += subt;
        }
      }
      fprintf(err, "\t%d\n", total);
    }
  }

  // print match/print counts
  printMatch(err, count, match, statusOpt);

  // print full breakdown
  if (err != NULL) {
    printSamples(err, elim, match, categ, opt, revmOpt);
    if (fclose(err))
      exit(error("", ERRCLOSE));
  }

  if (categ != NULL) {
    for (int i = 0; i < ERRCAT; i++)
      if (categ[i] != NULL)
        free(categ[i]);
    free(categ);
  }
}

/* void printFil()
 * Prints verbose filtering results.
 */
void printFil(FILE* fil, char* header, Midtag* mid,
    int outcome, char* categ, int flen, int tlen) {
  fprintf(fil, "%s\t%s\t%s\t", header, mid->name,
    mid->prim->name);
  if (outcome == TRUNC)
    fprintf(fil, "%s\t%s\t%d\t%d", TRUNCL, categ,
      flen, tlen);
  else if (outcome == ELIM)
    fprintf(fil, "%s\t%s\t%d", ELIML, categ, flen);
  else
    fprintf(fil, "%s\t\t\t%d", NEITHER, tlen);
  fprintf(fil, "\n");
}

/* void readFile()
 * Parses the sff.txt file.
 */
void readFile(FILE* file, FILE* fasta, int midMis,
    int primMis, int minSeqLength, int maxSeqLength,
    int maxTrunc, int maxAmbig, int okAmbig, int maxHomo,
    int okHomo, float avgQual, int windowLen, float windowAvg,
    int windowOpt, int minFlowLen, int maxFlows,
    float minNoise, float maxNoise, float maxValue,
    float noFlow, int seqCh, int qualCh, int flowCh,
    int samples, int midOpt, FILE* err, char** categ,
    float** inter, int denOpt, int denpOpt,
    FILE* denFasta, FILE* missFile, int** miss,
    int revmOpt, int revqOpt, int revMis, int trieOpt,
    int perOpt, float max, int statusOpt, FILE* fil) {

  // save header information
  int numReads = 0;
  if (statusOpt) {
    numReads = getInt(getData(file, line, NUMRE));
    printf("Filtering");
    if (denOpt)
      printf(" and Denoising (%s)", trieOpt ?
        "by trie" : "first iteration");
    printf(":\n");
  }
  int numFlows = getInt(getData(file, line, NUMFL));
  char* flowOrder = (char*) memalloc(numFlows + 1);
  char* orderLine = getData(file, line, CHARS);
  strcpy(flowOrder, orderLine);
  int reg = checkOrder(flowOrder, numFlows);
  for (Primer* p = primo; p != NULL; p = p->next) {
    fprintf(p->out, "%d %s\n", numFlows, flowOrder);
    if (denOpt && !trieOpt)
      p->prev = p->dummy;
  }
  numFlows++;

  // initialize counting variables
  int elim[samples][ETCAT][ERRCAT];
  int match[samples][ETCAT];
  for (int i = 0; i < samples; i++)
    for (int j = 0; j < ETCAT; j++) {
      match[i][j] = 0;
      for (int k = 0; k < ERRCAT; k++)
        elim[i][j][k] = 0;
    }
  int count = 0, print = 1;

  // initialize variables for saving data
  int maxlen = 2 * numFlows;
  char* header = (char*) memalloc(HEADER + 1);
  char* flowl = (char*) memalloc(MAX_SIZE);
  char* indexl = (char*) memalloc(MAX_SIZE);
  char* quall = (char*) memalloc(MAX_SIZE);
  float* flowgram = (float*) memalloc(numFlows * sizeof(float));
  int* index = (int*) memalloc(maxlen * sizeof(int));
  char* seq = (char*) memalloc(maxlen + 1);
  int* qual = (int*) memalloc(maxlen * sizeof(int));

  // read file
  while (fgets(line, MAX_SIZE, file) != NULL) {
    if (line[0] == '>') {
      count++;
      if (numReads)
        printf("%sCompleted %.1f%%\r", TAB,
          100.0f * count / numReads);

      // save header and sequence data
      int i;
      for (i = 0; line[i + 1] != '\n' && i < HEADER; i++)
        header[i] = line[i + 1];
      header[i] = '\0';
      int clipl = getInt(getData(file, line, CQL)) - 1;
      int clipr = getInt(getData(file, line, CQR));
      char* flows = getData(file, flowl, FLOWG);
      char* indexes = getData(file, indexl, FLOWI);
      char* seql = getData(file, line, BASE);

      // save bases -- after clipping
      int len = clipr - clipl;
      if (len > maxlen) {
        error(header, ERRLEN);
        len = maxlen;
      }
      for (i = 0; i < len; i++)
        seq[i] = seql[i + clipl + 1];
      seq[i] = '\0';

      // find matching mid tag - primer
      Midtag* mid = findMid(seq, midMis, primMis);
      if (mid != NULL) {

        match[mid->num][MATCH]++;
        int errNum = NOERR;
        int end = len;
        int min = strlen(mid->seq) + strlen(mid->prim->seq) + 1;  // default min. length
        int minLength = min > minSeqLength ? min : minSeqLength;

        // check sequence
        if (seqCh) {
          end = checkSeq(seq, len, minLength, maxSeqLength,
            maxTrunc, maxAmbig, okAmbig, maxHomo, okHomo,
            revmOpt, revqOpt, mid->prim->rev, revMis, &errNum);
          if (!end) {
            elim[mid->num][ELIM][errNum]++;
            if (fil != NULL)
              printFil(fil, header, mid, ELIM, categ[errNum],
                len, end);
            continue;
          }
        }

        // check quality scores
        if (qualCh) {
          char* quals = getData(file, quall, QUAL);
          loadInd(quals, qual, clipl, end, numFlows);

          int errNumQ = NOERR;
          int qualEnd = checkQual(qual, end, minLength,
            avgQual, windowLen, windowAvg, windowOpt, &errNumQ);
          if (!qualEnd) {
            elim[mid->num][ELIM][errNumQ]++;
            if (fil != NULL)
              printFil(fil, header, mid, ELIM, categ[errNumQ],
                len, end);
            continue;
          }
          if (qualEnd < end) {
            errNum = errNumQ;
            end = qualEnd;
          }
        }

        // save flowgram and flow indexes
        int begin = loadInd(indexes, index, clipl, end + 1, numFlows);
        int last = errNum && index[end - 1] < index[end] ?
          index[end] - 1 : index[end - 1]; // number of "good" flows
        int minFlows = index[minLength - 1] > minFlowLen ?
          index[minLength - 1] : minFlowLen; // min. flowgram length
        if (errNum && last < minFlows) {
          elim[mid->num][ELIM][errNum]++;
          if (fil != NULL)
            printFil(fil, header, mid, ELIM, categ[errNum],
              len, end);
          continue;
        }
        loadFlow(flows, flowgram, last, max);

        // check flowgram
        if (flowCh) {
          int errNumF = NOERR;
          int flowEnd = checkFlows(flowgram, begin, last, minFlows,
            maxFlows, minNoise, maxNoise, maxValue, noFlow, &errNumF);
          if (flowEnd < minFlows) {
            elim[mid->num][ELIM][errNumF]++;
            if (fil != NULL)
              printFil(fil, header, mid, ELIM, categ[errNumF],
                len, end);
            continue;
          }

          if (flowEnd < last) {
            errNum = errNumF;
            last = flowEnd;

            // determine new end
            int j;
            for (j = 0; index[j] <= flowEnd; j++) ;
            end = j;
          }
        }

        // print output
        match[mid->num][PRINT]++;
        elim[mid->num][TRUNC][errNum]++;
        if (fil != NULL)
          printFil(fil, header, mid, errNum ? TRUNC : -1,
            categ[errNum], len, end);
        int start = index[min - 2] - 1; // last flow value of primer
        seq[end] = '\0';
        printClean(fasta, header, mid, flowgram, start, last,
          seq, midOpt, reg, print);
        print++;

        // denoise flowgram
        if (denOpt) {
          if (trieOpt)
            trieStart(header, mid, flowgram, start, last,
              inter, reg, missFile, miss);
          else {
            Read* r = denoise(header, mid, flowgram, start,
              last, numFlows, inter, reg, missFile, miss);
            mid->prim->prev->next = r;
            mid->prim->prev = r;
          }
        }

      }
    }
  }

  // print cleaning statistics
  if (err != NULL || statusOpt)
    printStats(err, elim, match, categ, windowOpt,
      count, samples, revmOpt, statusOpt);

  // finish denoising (2nd iter [if nec.], print output)
  if (denOpt) {
    if (trieOpt)
      printDenTrie(denFasta, maxlen, flowOrder, denpOpt,
        midOpt, perOpt, reg, missFile, miss, max, statusOpt);
    else
      printDenoise(denFasta, flowOrder, denpOpt, midOpt,
        perOpt, reg, missFile, miss, inter, numFlows,
        max, statusOpt);

    // free inter array
    int maxm = (int) ((max + 0.01f) * 100.0f);
    for (int i = 0; i < maxm; i++)
      free(inter[i]);
    free(inter);
  }

  // free memory
  free(flowOrder);
  free(header);
  free(flowl);
  free(indexl);
  free(quall);
  free(seq);
  free(qual);
  free(flowgram);
  free(index);
  if (fclose(file) || (fasta != NULL && fclose(fasta))
      || (fil != NULL && fclose(fil)))
    exit(error("", ERRCLOSE));

}

/* char rc(char)
 * Returns the complement of the given base.
 */
char rc(char in) {
  char out;
  if (in == 'A')
    out = 'T';
  else if (in == 'T')
    out = 'A';
  else if (in == 'C')
    out = 'G';
  else if (in == 'G')
    out = 'C';
  else if (in == 'Y')
    out = 'R';
  else if (in == 'R')
    out = 'Y';
  else if (in == 'W')
    out = 'W';
  else if (in == 'S')
    out = 'S';
  else if (in == 'K')
    out = 'M';
  else if (in == 'M')
    out = 'K';
  else if (in == 'B')
    out = 'V';
  else if (in == 'V')
    out = 'B';
  else if (in == 'D')
    out = 'H';
  else if (in == 'H')
    out = 'D';
  else if (in == 'N')
    out = 'N';
  else
    exit(error("", ERRPRIM));
  return out;
}

/* char* revComp(char*)
 * Reverse-complements the given sequence.
 */
void revComp(char* out, char* seq) {
  int i = strlen(seq) - 1;
  int j;
  for (j = 0; i > -1; j++) {
    char nuc = seq[i--];
    out[j] = rc(nuc);
  }
  out[j] = '\0';
}

/* char* getMids()
 * Loads mid tags from the given file and saves
 *   them to the given primer.
 */
char* getMids(FILE* master, Primer* p, int revOpt,
    int* count) {
  Midtag* prev;
  while (fgets(line, MAX_SIZE, master) != NULL) {
    char* in = strtok(line, CSV);

    if (!strcmp(in, MIDTAG)) {
      char* name = strtok(NULL, CSV);
      char* seq = strtok(NULL, DELIM);
      if (name != NULL && seq != NULL) {

        // check for duplicate
        for (Midtag* mc = p->first; mc != NULL; mc = mc->next)
          if (!strcmp(mc->name, name))
            exit(error(name, ERRMREP));

        // create midtag
        Midtag* m = (Midtag*) memalloc(sizeof(Midtag));
        m->name = (char*) memalloc(1 + strlen(name));
        m->seq = (char*) memalloc(1 + strlen(seq));
        strcpy(m->name, name);
        strcpy(m->seq, seq);

        m->next = NULL;
        m->prim = p;
        m->num = (*count)++;

        if (p->first == NULL)
          p->first = m;
        else
          prev->next = m;
        prev = m;
      } else
        fprintf(stderr, "Error!  Skipping invalid line: %s,%s,%s\n",
          in, name, seq);

    } else if (!strcmp(in, REVERSE) && revOpt) {
      char* seq = strtok(NULL, DELIM);
      if (seq != NULL && p->rev == NULL) {
        p->rev = (char*) memalloc(1 + strlen(seq));
        revComp(p->rev, seq);
      } else
        fprintf(stderr, "Error!  Skipping invalid line: %s,%s\n",
          in, seq);

    } else if (!strcmp(in, PRIMER))
      return in;

  }
  return NULL;
}

/* FILE* overWrite(char*)
 * Prompts user to overwrite the given filename.
 */
FILE* overWrite(char* name) {
  error(name, ERREXIST);
  fprintf(stderr, "OK to overwrite? (y/n)  ");
  if (fgets(line, MAX_SIZE, stdin) != NULL &&
      (line[0] == 'y' || line[0] == 'Y'))
    return fopen(name, "w");
  return NULL;
}

/* FILE* openWrite(char*)
 * If the filename doesn't exist, opens it for writing.
 * Otherwise, prompts to overwrite.
 */
FILE* openWrite(char* outFile) {
  FILE* out = fopen(outFile, "r");
  if (out != NULL) {
    if (fclose(out))
      exit(error("", ERRCLOSE));
    out = overWrite(outFile);
  } else
    out = fopen(outFile, "w");
  if (out == NULL)
    exit(error(outFile, ERROPENW));
  return out;
}

/* char* getOutput(char*, char*)
 * Concatenates the 2nd parameter with the first.
 */
char* getOutput(char* in, char* ext) {
  char* out = (char*) memalloc(strlen(in) + strlen(ext) + 2);
  int i, j;
  for (i = 0; in[i] != '\0'; i++)
    out[i] = in[i];
  if (ext[0] != '.')
    out[i++] = '.';
  for (j = 0; ext[j] != '\0'; j++)
    out[i + j] = ext[j];
  out[i + j] = '\0';
  return out;
}

/* FILE* openPrimer(char*, char*, int)
 * Opens the primer's file for reading or writing
 *   depending on the opt (0 = "r", 1 = "w").
 */
FILE* openPrimer(char* name, char* ext, int opt) {
  char* file = getOutput(name, ext);
  FILE* out;
  if (opt)
    out = openWrite(file);
  else {
    out = fopen(file, "r");
    if (out == NULL)
      exit(error(file, ERROPENR));
  }
  free(file);
  return out;
}

/* char* findPrim(FILE*)
 * Finds the next primer line from the FILE*.
 */
char* findPrim(FILE* master) {
  while (fgets(line, MAX_SIZE, master) != NULL) {
    char* in = strtok(line, CSV);
    if (!strcmp(in, PRIMER))
      return in;
  }
  return NULL;
}

/* int loadSeqs(FILE*)
 * Loads the primers from the given file.
 * Returns the number of samples.
 */
int loadSeqs(FILE* master, char* flowExt, int opt,
    int denpOpt, char* denfExt, char* denmExt,
    int perOpt, char* perExt, char* chMapExt,
    int trieOpt, int revOpt) {

  char* in = findPrim(master);
  Primer* prev;
  int count = 0;
  while (in != NULL) {

    // load name and sequence
    char* name = strtok(NULL, CSV);
    char* seq = strtok(NULL, DELIM);
    if (name == NULL || seq == NULL) {
      fprintf(stderr, "Error!  Skipping invalid line: %s,%s,%s\n",
        in, name, seq);
      in = findPrim(master);
      continue;
    }

    // check for duplicate
    for (Primer* pc = primo; pc != NULL; pc = pc->next)
      if (!strcmp(pc->name, name))
        exit(error(name, ERRPREP));

    // create primer
    Primer* p = (Primer*) memalloc(sizeof(Primer));
    p->name = (char*) memalloc(1 + strlen(name));
    p->seq = (char*) memalloc(1 + strlen(seq));
    strcpy(p->name, name);
    strcpy(p->seq, seq);

    p->rev = NULL;
    p->first = NULL;
    p->next = NULL;
    p->head = NULL;
    p->dummy = (Read*) memalloc(sizeof(Read));
    p->out = openPrimer(name, flowExt, opt);

    if (trieOpt) {
      p->root = (Node*) memalloc(sizeof(Node));
      p->root->child = NULL;
      p->root->next = NULL;
      p->root->first = NULL;
      p->root->flow = NULL;
    } else
      p->root = NULL;

    if (denpOpt) {
      p->den = openPrimer(name, denfExt, denpOpt);
      p->map = trieOpt ? NULL :
        openPrimer(name, denmExt, denpOpt);
    } else
      p->den = p->map = NULL;

    if (perOpt) {
      p->per = openPrimer(name, perExt, perOpt);
      p->cmap = openPrimer(name, chMapExt, perOpt);
    } else
      p->per = p->cmap = NULL;

    if (primo == NULL)
      primo = p;
    else
      prev->next = p;
    prev = p;

    // load mid tags
    in = getMids(master, p, revOpt, &count);

    if (revOpt && p->rev == NULL)
      exit(error(p->name, ERRNOREV));
  }

  return count;
}

/* void printCat()
 * Saves the filtering criteria.
 */
void printCat(char** categ, int ecat, char* cat,
    void* val1, void* val2, int type) {
  categ[ecat] = (char*) memalloc(MIDPRIM);
  int i;
  for (i = 0; cat[i] != '\0'; i++)
    categ[ecat][i] = cat[i];
  if (val1 != NULL) {
    if (type == 1)
      sprintf(categ[ecat] + i, " (%d)", *(int*)val1);
    else if (type == 2)
      sprintf(categ[ecat] + i, " (%.2f)", *(float*)val1);
    else if (type == 3)
      sprintf(categ[ecat] + i, " (length=%d, qual=%.1f)",
        *(int*)val1, *(float*)val2);
    else if (type == 4)
      sprintf(categ[ecat] + i, " (%.2f-%.2f)",
        *(float*)val1, *(float*)val2);
  } else
    categ[ecat][i] = '\0';
}

/* char** getCateg()
 * Determines which filtering criteria are being used.
 */
char** getCateg(int minSeqLength, int maxSeqLength,
    int maxTrLength, int maxAmbig, int okAmbig,
    int maxHomo, int okHomo, float avgQual, int windowLen,
    float windowAvg, int minFlows, int maxFlows,
    float minNoise, float maxNoise, float maxValue,
    float noFlow, int revmOpt, int revqOpt) {

  char** categ = (char**) memalloc(ERRCAT * sizeof(char*));
  for (int i = 0; i < ERRCAT; i++)
    categ[i] = NULL;

  if (minSeqLength)
    printCat(categ, EMINSLEN, DMINSLEN, &minSeqLength, NULL, 1);
  if (maxSeqLength != MAX_SIZE)
    printCat(categ, EMAXSLEN, DMAXSLEN, &maxSeqLength, NULL, 1);
  if (maxTrLength != MAX_SIZE)
    printCat(categ, EMAXTRLEN, DMAXTRLEN, &maxTrLength, NULL, 1);
  if (maxAmbig != MAX_SIZE)
    printCat(categ, EMAXAMBIG, DMAXAMBIG, &maxAmbig, NULL, 1);
  if (okAmbig != MAX_SIZE)
    printCat(categ, EOKAMBIG, DOKAMBIG, &okAmbig, NULL, 1);
  if (maxHomo)
    printCat(categ, EMAXHOMO, DMAXHOMO, &maxHomo, NULL, 1);
  if (okHomo)
    printCat(categ, EOKHOMO, DOKHOMO, &okHomo, NULL, 1);
  if (revmOpt || revqOpt)
    printCat(categ, EREVERSE, DREVERSE, NULL, NULL, 0);
  if (avgQual)
    printCat(categ, EAVGQUAL, DAVGQUAL, &avgQual, NULL, 2);
  if (windowAvg)
    printCat(categ, EWINDOW, DWINDOW, &windowLen, &windowAvg, 3);
  if (minFlows)
    printCat(categ, EMINFLEN, DMINFLEN, &minFlows, NULL, 1);
  if (maxFlows != MAX_SIZE)
    printCat(categ, EMAXFLEN, DMAXFLEN, &maxFlows, NULL, 1);
  if (minNoise)
    printCat(categ, EFLOWINT, DFLOWINT, &minNoise, &maxNoise, 4);
  if (maxValue != DEFMAXFLOW)
    printCat(categ, EMAXNVAL, DMAXNVAL, &maxValue, NULL, 2);
  if (noFlow)
    printCat(categ, ENOFLOW, DNOFLOW, &noFlow, NULL, 2);

  return categ;
}

/* void loadDist()
 * Loads asymmetric distances from given file.
 */
void loadDist(FILE* sd, float** inter, float val,
    int max) {
  for (int i = 0; i < max; i++) {
    float f = i / 100.0f;

    // load values from stddev file
    do {
      if (fgets(line, MAX_SIZE, sd) == NULL)
        exit(error("", ERRSDVAL));
    } while (line[0] == '#' || line[0] == '\n');
    char* dist1 = strtok(line, CSV);
    char* dist2 = strtok(NULL, DELIM);
    if (dist1 == NULL || dist2 == NULL)
      exit(error("", ERRSDVAL));

    float d1 = val * getFloat(dist1);
    float d2 = val * getFloat(dist2);
    inter[i][0] = f - d1;
    inter[i][1] = f + d2;
  }
}

/* int checkSD()
 * Checks file of distances to determine if they
 *   are asymmetric.
 */
int checkSD(FILE* sd, float** inter, float val,
    int max) {
  if (sd != NULL) {
    do {
      if (fgets(line, MAX_SIZE, sd) == NULL)
        exit(error("", ERRSDVAL));
    } while (line[0] == '#' || line[0] == '\n');
    for (int i = 0; line[i] != '\0'; i++)
      if (line[i] == ',' || line[i] == '\t') {
        rewind(sd);
        loadDist(sd, inter, val, max);
        return 0;
      }
    rewind(sd);
  }
  return 1;
}

/* void makeInter(FILE*, float[][], float)
 * Makes the flow intervals for each flow value.
 */
void makeInter(FILE* sd, float** inter, float val,
    int max) {

  if (checkSD(sd, inter, val, max)) {
    for (int i = 0; i < max; i++) {
      float f = i / 100.0f;
      float diff = val;

      // load value from stddev file
      if (sd != NULL) {
        do {
          if (fgets(line, MAX_SIZE, sd) == NULL)
            exit(error("", ERRSDVAL));
        } while (line[0] == '#' || line[0] == '\n');
        float fl = getFloat(line);
        diff *= fl;
      }

      inter[i][0] = f - diff;
      inter[i][1] = f + diff;
    }
  }
}

/* int defDenoise(float, float)
 * Returns denOpt if no denoising parameter is specified
 *   and default (cleaning and denoising) is used.
 */
int defDenoise(float cons, float perc) {
  if (!cons && !perc) {
    fprintf(stderr, "Warning!  No denoising parameter specified\n");
    fprintf(stderr, "Do you want to denoise? (y/n)  ");
    if (fgets(line, MAX_SIZE, stdin) != NULL &&
        (line[0] == 'y' || line[0] == 'Y'))
      exit(error("", ERRDEN));
    return 0;
  }
  return 1;
}

/* void checkParams()
 * Check parameters for issues.
 */
void checkParams(char* masterFile, FILE** sff, char* sffFile,
    FILE** out, char* outFile, char* flowExt, FILE** err,
    char* errFile, char* sdFile, int denpOpt, char* denfExt,
    char* denmExt, FILE** fasta, char* fastaFile, FILE** misses,
    char* missFile, FILE** fil, char* filFile,
    char*** categ, int midMis, int primMis,
    int minSeqLength, int maxSeqLength, int maxTrLength,
    int maxAmbig, int okAmbig, int maxHomo, int okHomo,
    float avgQual, int windowLen, float windowAvg,
    int minFlows, int maxFlows, float minNoise,
    float maxNoise, float maxValue, float noFlow,
    int* seqCh, int* qualCh, int* flowCh, int* count,
    float*** inter, int* cleanOpt, int* denOpt,
    float consInt, float percInt, int chimOpt, char* chimExt,
    char* chMapExt, int revmOpt, int revqOpt, int revMis,
    int trieOpt, float maxFlow, int*** miss) {

  // open master file
  if (masterFile == NULL)
    exit(error("master.csv", ERRNEED));
  FILE* master = fopen(masterFile, "r");
  if (master == NULL)
    exit(error(masterFile, ERROPENR));

  if (!*cleanOpt && !*denOpt) {
    *cleanOpt = 1;
    *denOpt = defDenoise(consInt, percInt);
  }

  // check cleaning options
  if (*cleanOpt) {

    // open sff.txt file
    if (sffFile == NULL)
      exit(error("sff.txt", ERRNEED));
    *sff = fopen(sffFile, "r");
    if (*sff == NULL)
      exit(error(sffFile, ERROPENR));

    // sequences
    if (revmOpt && revqOpt)
      exit(error("", ERRREV));
    if (midMis < 0 || primMis < 0 ||
        ((revmOpt || revqOpt) && revMis < 0))
      exit(error("", ERRMISM));
    if (minSeqLength < 0 || minSeqLength > maxSeqLength
        || minSeqLength > maxTrLength)
      exit(error("", ERRSLEN));
    if (!minSeqLength && maxSeqLength == MAX_SIZE &&
        maxTrLength == MAX_SIZE && maxAmbig == MAX_SIZE &&
        okAmbig == MAX_SIZE && !maxHomo && !okHomo &&
        !revmOpt && !revqOpt)
      *seqCh = 0;

    // quality scores
    if ((windowLen && !windowAvg) || (!windowLen && windowAvg)
        || windowLen < 0)
      exit(error("", ERRWIND));
    if (!avgQual && !windowAvg)
      *qualCh = 0;

    // flowgrams
    if (minNoise < 0.0f || maxNoise > maxFlow || minNoise > maxNoise
        || (minNoise && !maxNoise) || (!minNoise && maxNoise))
      exit(error("", ERRINTER));
    if (minFlows < 0 || minFlows > maxFlows)
      exit(error("", ERRFLEN));
    if (maxValue > DEFMAXFLOW)
      exit(error("", ERRMAXFL));
    if (!minNoise && maxValue == DEFMAXFLOW &&
        !minFlows && maxFlows == MAX_SIZE && !noFlow)
      *flowCh = 0;

    // open output files
    if (outFile != NULL)
      *out = openWrite(outFile);
    if (filFile != NULL || errFile != NULL) {
      if (filFile != NULL) {
        *fil = openWrite(filFile);
        fprintf(*fil, "%s", FILHEAD);
      }
      if (errFile != NULL)
        *err = openWrite(errFile);
      *categ = getCateg(minSeqLength, maxSeqLength,
        maxTrLength, maxAmbig, okAmbig, maxHomo, okHomo,
        avgQual, windowLen, windowAvg, minFlows, maxFlows,
        minNoise, maxNoise, maxValue, noFlow, revmOpt,
        revqOpt);
    }
  }

  // set default parameters
  if (flowExt == NULL)
    flowExt = DEFFLOWEXT;
  if (*denOpt) {
    if (denpOpt) {
      if (denfExt == NULL)
        denfExt = DEFFEXT;
      if (denmExt == NULL)
        denmExt = DEFMEXT;
    }
    if (chimOpt && chimExt == NULL)
      chimExt = DEFCHEXT;
    if (chimOpt && chMapExt == NULL)
      chMapExt = DEFCHMAP;
  } else
    denpOpt = chimOpt = 0;

  // load primer and mid tags from master file
  *count = loadSeqs(master, flowExt, *cleanOpt,
    denpOpt, denfExt, denmExt, chimOpt, chimExt,
    chMapExt, trieOpt, revmOpt || revqOpt);
  if (fclose(master))
    exit(error("", ERRCLOSE));

  // check denoising options
  if (*denOpt) {

    if ((consInt && percInt) || (!consInt && !percInt))
      exit(error("", ERRDEN));
    if (consInt < 0.0f || consInt > DEFMAXFLOW || percInt < 0.0f)
      exit(error("", ERRDENVAL));

    if (fastaFile == NULL)
      fastaFile = DEFDENFILE;
    *fasta = openWrite(fastaFile);

    // check maxFlow
    if (maxFlow < 0)
      exit(error("", ERRMAXF));
    int max = (int) ((maxFlow + 0.01f) * 100.0f);

    // make interval array
    *inter = (float**) memalloc(max * sizeof(float*));
    for (int i = 0; i < max; i++)
      (*inter)[i] = (float*) memalloc(ETCAT * sizeof(float));
    if (percInt) {
      if (sdFile == NULL)
        sdFile = DEFSDFILE;
      FILE* sd = fopen(sdFile, "r");
      if (sd == NULL)
        exit(error(sdFile, ERROPENR));
      makeInter(sd, *inter, percInt, max);
      if (fclose(sd))
        exit(error("", ERRCLOSE));
    } else
      makeInter(NULL, *inter, consInt, max);

    // make "misses" array
    if (missFile != NULL) {
      *misses = openWrite(missFile);
      *miss = (int**) memalloc(max * sizeof(int*));
      for (int i = 0; i < max; i++) {
        (*miss)[i] = (int*) memalloc(max * sizeof(int));
        for (int j = 0; j < max; j++)
          (*miss)[i][j] = 0;
      }
    }

  }

}

/* void getParams(int, char**)
 * Gets the parameters from the command line.
 */
void getParams(int argc, char* argv[]) {

  // initialize parameters
  char* master = NULL, *sffFile = NULL, *outFile = NULL,
    *flowExt = NULL, *errFile = NULL, *sdFile = NULL,
    *denfExt = NULL, *denmExt = NULL, *fastaFile = NULL,
    *missFile = NULL, *chimExt = NULL, *chMapExt = NULL,
    *filFile = NULL;
  int midOpt = 0, windowOpt = 0, cleanOpt = 0, denOpt = 0,
    denpOpt = 0, chimOpt = 0, perOpt = 0,
    revmOpt = 0, revqOpt = 0, trieOpt = 0,
    statusOpt = 0;
  int midMis = 0, primMis = 0, minSeqLength = 0,
    maxSeqLength = MAX_SIZE, maxTrLength = MAX_SIZE,
    maxAmbig = MAX_SIZE, okAmbig = MAX_SIZE, maxHomo = 0,
    okHomo = 0, windowLen = 0, minFlows = 0,
    maxFlows = MAX_SIZE, revMis = 0;
  float avgQual = 0.0f, windowAvg = 0.0f, minNoise = 0.0f,
    maxNoise = 0.0f, maxValue = DEFMAXFLOW, noFlow = 0.0f,
    consInt = 0.0f, percInt = 0.0f, maxFlow = DEFMAXFLOW;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (!strcmp(argv[i], WINDOWOPT))
      windowOpt = 1;
    else if (!strcmp(argv[i], NOMIDOPT))
      midOpt = 1;
    else if (!strcmp(argv[i], CLEANOPT))
      cleanOpt = 1;
    else if (!strcmp(argv[i], DENOPT))
      denOpt = 1;
    else if (!strcmp(argv[i], BOTHOPT))
      cleanOpt = denOpt = 1;
    else if (!strcmp(argv[i], DENPOPT))
      denpOpt = 1;
    else if (!strcmp(argv[i], CHIMOPT))
      chimOpt = 1;
    else if (!strcmp(argv[i], REVMOPT))
      revmOpt = 1;
    else if (!strcmp(argv[i], REVQOPT))
      revqOpt = 1;
    else if (!strcmp(argv[i], TRIEOPT))
      trieOpt = 1;
    else if (!strcmp(argv[i], STATUSOPT))
      statusOpt = 1;
    else if (i < argc - 1) {
      if (!strcmp(argv[i], MASTERFILE))
        master = argv[++i];
      else if (!strcmp(argv[i], SFFFILE))
        sffFile = argv[++i];
      else if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], FLOWEXT))
        flowExt = argv[++i];
      else if (!strcmp(argv[i], ERRFILE))
        errFile = argv[++i];
      else if (!strcmp(argv[i], SDFILE))
        sdFile = argv[++i];
      else if (!strcmp(argv[i], MISSFILE))
        missFile = argv[++i];
      else if (!strcmp(argv[i], DENFEXT))
        denfExt = argv[++i];
      else if (!strcmp(argv[i], DENMEXT))
        denmExt = argv[++i];
      else if (!strcmp(argv[i], DENFASTA))
        fastaFile = argv[++i];
      else if (!strcmp(argv[i], FILFILE))
        filFile = argv[++i];
      else if (!strcmp(argv[i], PEREXT)) {
        chimExt = argv[++i];
        perOpt = 1;
      } else if (!strcmp(argv[i], UCHEXT)) {
        chimExt = argv[++i];
        perOpt = 0;
      } else if (!strcmp(argv[i], UMAPEXT))
        chMapExt = argv[++i];
      else if (!strcmp(argv[i], MIDMIS))
        midMis = getInt(argv[++i]);
      else if (!strcmp(argv[i], PRIMMIS))
        primMis = getInt(argv[++i]);
      else if (!strcmp(argv[i], MINSLEN))
        minSeqLength = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXSLEN))
        maxSeqLength = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXTRLEN))
        maxTrLength = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXAMBIG))
        maxAmbig = getInt(argv[++i]);
      else if (!strcmp(argv[i], OKAMBIG))
        okAmbig = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXHOMO))
        maxHomo = getInt(argv[++i]);
      else if (!strcmp(argv[i], OKHOMO))
        okHomo = getInt(argv[++i]);
      else if (!strcmp(argv[i], AVGQUAL))
        avgQual = getFloat(argv[++i]);
      else if (!strcmp(argv[i], WINDOWLEN))
        windowLen = getInt(argv[++i]);
      else if (!strcmp(argv[i], WINDOWAVG))
        windowAvg = getFloat(argv[++i]);
      else if (!strcmp(argv[i], MINFLEN))
        minFlows = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXFLEN))
        maxFlows = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXFLOW))
        maxFlow = getFloat(argv[++i]);
      else if (!strcmp(argv[i], MININT))
        minNoise = getFloat(argv[++i]);
      else if (!strcmp(argv[i], MAXINT))
        maxNoise = getFloat(argv[++i]);
      else if (!strcmp(argv[i], MAXNVAL))
        maxValue = getFloat(argv[++i]);
      else if (!strcmp(argv[i], NOFLOW))
        noFlow = getFloat(argv[++i]);
      else if (!strcmp(argv[i], CINTER))
        consInt = getFloat(argv[++i]);
      else if (!strcmp(argv[i], ZINTER))
        percInt = getFloat(argv[++i]);
      else if (!strcmp(argv[i], REVMIS))
        revMis = getInt(argv[++i]);
      else
        exit(error(argv[i], ERRINVAL));
    } else
      exit(error(argv[i], ERRINVAL));
  }

  // check parameters
  FILE* sff = NULL, *out = NULL, *err = NULL, *fasta = NULL,
    *misses = NULL, *fil = NULL;
  int seqCh = 1, qualCh = 1, flowCh = 1; // booleans for checking
  int samples = 0; // number of samples analyzed
  char** categ = NULL; // categories of elim/trunc criteria
  float** inter;
  int** miss;

  checkParams(master, &sff, sffFile, &out, outFile,
    flowExt, &err, errFile, sdFile, denpOpt, denfExt,
    denmExt, &fasta, fastaFile, &misses, missFile,
    &fil, filFile, &categ, midMis, primMis,
    minSeqLength, maxSeqLength, maxTrLength,
    maxAmbig, okAmbig, maxHomo, okHomo, avgQual,
    windowLen, windowAvg, minFlows, maxFlows,
    minNoise, maxNoise, maxValue, noFlow,
    &seqCh, &qualCh, &flowCh,
    &samples, &inter, &cleanOpt, &denOpt, consInt, percInt,
    chimOpt, chimExt, chMapExt, revmOpt, revqOpt, revMis,
    trieOpt, maxFlow, &miss);

  // initialize cleaning/denoising
  if (cleanOpt)
    readFile(sff, out, midMis, primMis, minSeqLength,
      maxSeqLength, maxTrLength, maxAmbig, okAmbig, maxHomo,
      okHomo, avgQual, windowLen, windowAvg, windowOpt,
      minFlows, maxFlows, minNoise, maxNoise, maxValue,
      noFlow, seqCh, qualCh, flowCh, samples, midOpt, err,
      categ, inter, denOpt, denpOpt, fasta, misses, miss,
      revmOpt, revqOpt, revMis, trieOpt, perOpt, maxFlow,
      statusOpt, fil);

  else
    readClean(midOpt, inter, denpOpt, fasta, misses, miss,
      trieOpt, perOpt, maxFlow, statusOpt);

}

/* int main(int, char*)
 * Main.
 */
int main(int argc, char* argv[]) {

  line = (char*) memalloc(MAX_SIZE);
  primo = NULL;

  getParams(argc, argv);
  freeMemory();

  return 0;
}
