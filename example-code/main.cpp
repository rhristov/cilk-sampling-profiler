// Copyright 2014

#include <assert.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <algorithm>
#include <cmath>
#include <ctime>

#include <sstream>
#include <string>
#include <vector>

#include "./tfk_utils.h"
//#include "sprp64_sf.h"
#include "sprp64.h"
#include "sprp32.h"
#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <cilk/reducer_opadd.h>
#include <google/profiler.h>
/*
std::string zToString(const ZZ &z) {
    std::stringstream buffer;
    buffer << z;
    return buffer.str();
}
*/
// Static variables for convenience.
// Most of these aren't needed, but are good for sanity checking.
static int total_updates = 0;
static int* starts;
static int* ends;
static int* success_vector;


uint64_t modular_pow(uint64_t base, uint64_t exponent, int modulus)
{
    uint64_t result = 1;
    while (exponent > 0)
    {
        if (exponent % 2 == 1)
            result = (result * base) % modulus;
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }
    return result;
}
int mul_inv(int a, int b)
{
	int b0 = b, t, q;
	int x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

int legendre(int x,int p)
{ // finds (x/p) as 1 or -1
        int m,k,p8,t;
 
        m=0;
        while(p>1)
        { /* main loop */
 
// extract powers of 2
                for (k=0;x%2==0;k++) x/=2;
                p8=p%8;
                if (k%2==1) m+=(p8*p8-1)/8;
 
// quadratic reciprocity
                t=p;  t%=x;
                p=x;  x=t;
                m+=(p8-1)*(x%4-1)/4;
                m%=2;
        }
        if (m==0) return 1;
        else      return (-1);
}
/* 
void find_primes(int start, int end) {

  if (start%2==0) start++;

  int* composites = (int*) calloc(64, sizeof(int));
  //int* composite_buffer = (int*) calloc(64, sizeof(int));
  //int* composite_buffer_last_n = (int*) calloc(64, sizeof(int));

  srand(0);
  int num_composite_found = 0;
  while(num_composite_found < 64) {
    uint32_t c = (uint32_t)rand();
    c = c%start;
    if (c % 2 != 0) {
      composites[num_composite_found] = c;
      //composite_buffer[num_composite_found] = pow(c, (start-1)/2);
      //composite_buffer_last_n[num_composite_found] = pow(c, (start-1)/2);
      num_composite_found++;
      printf("composite found %d\n", c);
    }
  }

  printf("trying to find primes now\n");
  for (int p = start; p < end; p++) {
    if (p % 2 == 0) continue;
    bool prime = true;
    // check phase
    for (int i = 0; i < 64; i++) {
      uint64_t result = modular_pow(composites[i], (p-1)/2, p);
      uint64_t result2 = ZZ::Jacobi(composites[i], p) % p;
      //printf("result %d, result2 %d\n", result, result2);
      if (result == result2) continue;
      prime = false;
      break;
    }
    if (prime) printf("found a prime %d\n", p);
  }

}
*/

const uint64_t bases[] = {2,325,9375,28178,450775,9780504,1795265022};

const uint32_t bases32[] = {15, 176006322, 4221622697};
const int bases32_cnt = 3;

const uint32_t bases321[] = {921211727};
const uint32_t bases321_limit = 49141;
const int bases321_cnt = 1;

const uint32_t bases322[] = {1143370, 2350307676};
const uint32_t bases322_limit = 360018361;
const int bases322_cnt = 2;

/*
const uint32_t bases323[] = {15, 176006322, 4221622697};
const uint64_t bases323_limit = ;
const int bases323_cnt = 3;
*/

const uint64_t bases1[] = {9345883071009581737};
const uint64_t bases1_limit = 341531;
const int bases1_cnt = 1;

const uint64_t bases2[] = {336781006125, 9639812373923155};
const uint64_t bases2_limit = 1050535501;
const int bases2_cnt = 2;

const uint64_t bases3[] = {4230279247111683200, 14694767155120705706, 16641139526367750375};
const uint64_t bases3_limit = 350269456337;
const int bases3_cnt = 3;

//const int bases_cnt = 1;

//const uint64_t bases[1] = {2};
const int bases_cnt = 7;

const uint64_t small_primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
const int small_primes_cnt = 19;

static char* sieve;

inline bool comp_filter(uint64_t i) {
//  return sieve[i] == 1;

  if (i%3==0) return true;
  if (i%5==0) return true;
  if (i%7==0) return true;
  if (i%11==0) return true;
  if (i%13==0) return true;
  if (i%17==0) return true;
  if (i%19==0) return true;
  if (i%23==0) return true;
  if (i%29==0) return true;
  if (i%31==0) return true;
  if (i%37==0) return true;
  if (i%41==0) return true;
  if (i%43==0) return true;
  if (i%47==0) return true;
  if (i%53==0) return true;
  if (i%59==0) return true;
  if (i%61==0) return true;
  if (i%67==0) return true;
// NOTE(TFK): Tests below aren't worth it.
/*
  if (i%71==0) return true;
  if (i%73==0) return true;
  if (i%79==0) return true;
  if (i%83==0) return true;
  if (i%89==0) return true;
  if (i%97==0) return true;
  if (i%101==0) return true;
  if (i%103==0) return true;
  if (i%107==0) return true;
  if (i%109==0) return true;
  if (i%113==0) return true;
  if (i%127==0) return true;
  if (i%131==0) return true;
  if (i%137==0) return true;
  if (i%139==0) return true;
  if (i%149==0) return true;
  if (i%151==0) return true;
  if (i%157==0) return true;
  if (i%163==0) return true;
  if (i%167==0) return true;
  if (i%173==0) return true;
  if (i%179==0) return true;
  if (i%181==0) return true;
  if (i%191==0) return true;
  if (i%193==0) return true;
  if (i%197==0) return true;
  if (i%199==0) return true;
*/

  return false;
}


char* prime_s(uint64_t _start, uint64_t _end) {
  char* s = (char*) calloc(_end-_start + 1, sizeof(char));

  uint64_t primes[] = {3,5};
  uint64_t primes_cnt = 2;

  uint64_t index = _start;
  for (int i = _start; i < _end; i++) {
      if (s[i] == 0) {
        for (int j = 0; j < primes_cnt; j++) {
          if (i%primes[j]==0) {
            for (int k = _start; k < _end; k+=primes[j]) {
              s[k] = 1;
            }
          }
        }
      }
  }
  return s;
}
#define PARALLEL

void find_primes32(uint32_t _start, uint32_t _end) {
  /*ZZ start;
  start = _start;
  ZZ end;
  end = _end;*/
  uint32_t start = _start;
  uint32_t end = _end;
  uint32_t prime_count = 0;
  uint32_t invoke_count = 0;

  if (start%2==0) start++;
  if (start < bases321_limit) {
    uint32_t limit = bases321_limit;
    if (limit > end) limit = end;
    for (uint32_t i = start; i < limit; i += 2) {
      //if (i%3==0) continue;
      //if (i%5==0) continue;
      const uint32_t p = i;
      if (efficient_mr32(bases321, bases321_cnt, p)) {
        #ifdef PARALLEL
        __sync_fetch_and_add(&prime_count, 1);
        #else
        prime_count++;
        #endif
      }
        invoke_count++;
      //if (i%3==1) i+=2;
    }
    start = bases321_limit;
  }
  if (start%2==0) start++;
  if (start < bases322_limit) {
    uint64_t limit = bases322_limit;
    if (limit > end) limit = end;
    for (uint64_t i = start; i < limit; i += 2) {
      if (comp_filter(i)) continue;
      const uint64_t p = i;
      if (efficient_mr32(bases322, bases322_cnt, p)) {
        #ifdef PARALLEL
        __sync_fetch_and_add(&prime_count, 1);
        #else
        prime_count++;
        #endif
      }
        invoke_count++;
      //if (i%3==1) i+=2;
    }
    start = bases322_limit;
  }
/*
  if (start%2==0) start++;
  if (start < bases323_limit) {
    uint64_t limit = bases323_limit;
    if (limit > end) limit = end;
    for (uint64_t i = start; i < limit; i += 2) {
      if (comp_filter(i)) continue;
      const uint64_t p = i;
      if (efficient_mr64(bases323, bases323_cnt, p)) {
        #ifdef PARALLEL
        __sync_fetch_and_add(&prime_count, 1);
        #else
        prime_count++;
        #endif
        //__sync_fetch_and_add(&prime_count, 1);
      }
    }
    start = bases323_limit;
  }
*/
  if (start%2==0) start++;
  if (start < end) {
  for(uint32_t i = start; i < end; i += 2) {
    if (comp_filter(i)) continue;
    const uint32_t p = i;
    if (efficient_mr32(bases32, bases32_cnt, p)) {
        #ifdef PARALLEL
        __sync_fetch_and_add(&prime_count, 1);
        #else
        prime_count++;
        #endif
      //__sync_fetch_and_add(&prime_count, 1);
    }
        invoke_count++;
      //if (i%3==1) i+=2;
  }
    printf("doing slow one\n");
  }
/*
  while (start < end) {
    if (efficient_mr64(bases, bases_cnt, start)) {
      prime_count++;
    }
    start++;
    //NextPrime(start, start+1, 1);
    //prime_count++;
  }*/
  //prime_count--;
  printf("prime count is %llu\n", prime_count);
  printf("R-M invoke count is %llu\n", invoke_count);
}

#define PARALLEL
void find_primes(uint64_t _start, uint64_t _end) {
  /*ZZ start;
  start = _start;
  ZZ end;
  end = _end;*/
  uint64_t start = _start;
  uint64_t end = _end;
  #ifdef PARALLEL
  cilk::reducer_opadd<uint64_t> prime_count(0);
  #else
  uint64_t prime_count = 0;
  #endif

  if (start%2==0) start++;
  if (start < bases1_limit) {
    uint64_t limit = bases1_limit;
    if (limit > end) limit = end;
    cilk_for (uint64_t i = start; i < limit; i += 2) {
      //if (i%3==0) continue;
      //if (i%5==0) continue;
      const uint64_t p = i;
      if (efficient_mr64(bases1, bases1_cnt, p)) {
        #ifdef PARALLEL
        //__sync_fetch_and_add(&prime_count, 1);
        *prime_count += 1;
        #else
        prime_count++;
        #endif
      }
    }
    start = bases1_limit;
  }
  if (start%2==0) start++;
  if (start < bases2_limit) {
    uint64_t limit = bases2_limit;
    if (limit > end) limit = end;
    cilk_for (uint64_t i = start; i < limit; i += 2) {
      if (comp_filter(i)) continue;
      const uint64_t p = i;
      if (efficient_mr64(bases2, bases2_cnt, p)) {
        #ifdef PARALLEL
        //__sync_fetch_and_add(&prime_count, 1);
        *prime_count += 1;
        #else
        prime_count++;
        #endif
      }
    }
    start = bases2_limit;
  }
  if (start%2==0) start++;
  if (start < bases3_limit) {
    uint64_t limit = bases3_limit;
    if (limit > end) limit = end;
    cilk_for (uint64_t i = start; i < limit; i += 2) {
      if (comp_filter(i)) continue;
      const uint64_t p = i;
      if (efficient_mr64(bases3, bases3_cnt, p)) {
        #ifdef PARALLEL
        //__sync_fetch_and_add(&prime_count, 1);
        *prime_count += 1;
        #else
        prime_count++;
        #endif
        //__sync_fetch_and_add(&prime_count, 1);
      }
    }
    start = bases3_limit;
  }

  if (start%2==0) start++;
  if (start < end) {
  cilk_for(uint64_t i = start; i < end; i += 2) {
    if (comp_filter(i)) continue;
    const uint64_t p = i;
    if (efficient_mr64(bases, bases_cnt, p)) {
        #ifdef PARALLEL
        //__sync_fetch_and_add(&prime_count, 1);
        *prime_count += 1;
        #else
        prime_count++;
        #endif
      //__sync_fetch_and_add(&prime_count, 1);
    }
  }
    printf("doing slow one\n");
  }
/*
  while (start < end) {
    if (efficient_mr64(bases, bases_cnt, start)) {
      prime_count++;
    }
    start++;
    //NextPrime(start, start+1, 1);
    //prime_count++;
  }*/
  //prime_count--;
  printf("prime count is %llu\n", prime_count.get_value());
}



int main(int argc, char **argv) {
  uint64_t N = 1000;
  uint64_t REPS = 2000;
  bool VERIFY = false;
  if (argc > 2) {
    printf("N=%s REPS=%s \n", argv[1], argv[2]);
    N = StringToNumber<uint64_t>(std::string(argv[1]));
    REPS = StringToNumber<uint64_t>(std::string(argv[2]));
    if (argc > 3) {
      printf("Running affinity loop.\n");
      VERIFY = true;
    } else {
      printf("Running regular for loop\n");
    }
  } else {
    printf("Usage: ./main <N> <REPS> <?> ---"
        "include optional 3rd argument to run affinity loop\n");
    return 0;
  }
  //sieve = prime_s(N, REPS);
  //printf("done computing sieve\n");
  //char* x = (char*) malloc(sizeof(char) * (1<<40));
  //printf("x 0 is %d\n", x[0]);
/*  if (REPS < 4294967295) {
    find_primes32((uint32_t)N, (uint32_t)REPS);
  } else {
    find_primes32((uint32_t)N, 4294967295);
    if (N > 4294967295) {
      find_primes(4294967295, REPS);
    }
  }
  */
  ProfilerStart("profile.data"); 
      find_primes(N, REPS);
  ProfilerStop();
  return 0;
}
