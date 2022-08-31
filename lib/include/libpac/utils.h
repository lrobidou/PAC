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

uint64_t str2num(const std::string& str);

uint64_t revhash(uint64_t x);

uint64_t unrevhash(uint64_t x);

// TODO CAN BE IMPROVED
// rcb: reverse complement, binary
uint64_t rcb(uint64_t min, uint64_t n);

uint64_t canonize(uint64_t x, uint64_t n);

void createPathIfNotExists(const std::string& path);

uint64_t nuc2int(char c);

std::string intToString(uint64_t n);
uint64_t asm_log2(const uint64_t x);
uint64_t approx_power2(uint64_t n);

uint64_t nuc2intrc(char c);

int directory_exists(std::string& path);

uint64_t str2numstrand(const std::string& str);

uint64_t xorshift64(uint64_t& state);

uint64_t hash64shift(uint64_t key);

// max is excluded
uint64_t hash_family(uint64_t key, uint b, uint64_t max);

void Biogetline(zstr::ifstream* in, std::string& result, char type, uint K);

std::tuple<std::string, std::string> Biogetline(zstr::ifstream* in, char type, uint K);

char get_data_type(const std::string& filename);

bool exists_test(const std::string& name);

uint64_t getMemorySelfMaxUsed();

#endif
