#pragma once
#ifndef RAND_H
#define RAND_H

#include <math.h>
#include <ctime>
#include <random> 
#include "stdc++.h"
#include <algorithm>
using namespace std;

double _rand(double min, double max) { //[min,max)
    static thread_local random_device rd;
    static thread_local mt19937_64 generator(rd());
    uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}
vector<int> random_sample(int idx, int N, int g) {
    int k = g - 1;
    int population_size = N - idx;
    vector<int> population(population_size);
    iota(population.begin(), population.end(), idx);
    static random_device rd;
    static mt19937 gen(rd());
    shuffle(population.begin(), population.end(), gen);
    population.resize(k);
    return population;
}
vector<int> make_sorted_edge(int idx, int N, int g) {
    vector<int> sampled = random_sample(idx + 1, N, g);
    sampled.push_back(idx);
    sort(sampled.begin(), sampled.end());
    return sampled;
}
void shuffle(std::vector<int>& vec)
{
    // 高质量随机数引擎（Mersenne Twister）
    static thread_local std::mt19937 rng(
        std::random_device{}()
    );

    // Fisher–Yates 洗牌（标准实现）
    for (int i = static_cast<int>(vec.size()) - 1; i > 0; --i)
    {
        std::uniform_int_distribution<int> dist(0, i);
        std::swap(vec[i], vec[dist(rng)]);
    }
}
#endif