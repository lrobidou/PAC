#ifndef UTILH
#define UTILH

#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <zstr.hpp>

using namespace std;

uint64_t nuc2int(char c);

string intToString(uint64_t n);
uint64_t asm_log2(const uint64_t x);
uint64_t approx_power2(uint64_t n);

uint64_t nuc2intrc(char c);

int directory_exists(string& path);

uint64_t str2numstrand(const string& str);

uint64_t xorshift64(uint64_t& state);

uint64_t hash64shift(uint64_t key);

uint64_t hash_family(uint64_t key, uint b);

void Biogetline(zstr::ifstream* in, string& result, char type, uint K);

void Biogetline(zstr::ifstream* in, string& result, char type, uint K, string& header);

char get_data_type(const string& filename);

bool exists_test(const std::string& name);

uint64_t getMemorySelfMaxUsed();

#endif
