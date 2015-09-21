#pragma once
#include<vector>
#include<map>
#include<string>
#include<set>
#include<math.h>
#include<algorithm>
#include<fstream>
#include<iostream>
#include "classObject.h"
using namespace std;
typedef pair<int, double> PAIR;
//字符串分割函数
std::vector<std::string> split(std::string str, std::string pattern)
{
	std::string::size_type pos;
	std::vector<std::string> result;
	str += pattern;//扩展字符串以方便操作
	int size = str.size();

	for (int i = 0; i<size; i++)
	{
		pos = str.find(pattern, i);
		if (pos<size)
		{
			std::string s = str.substr(i, pos - i);
			result.push_back(s);
			i = pos + pattern.size() - 1;
		}
	}
	return result;
}
vector<vector<int>> statisticRate(string filename)
{
	ifstream userInfoFilePtr;
	userInfoFilePtr.open(filename);
	string temp;
	vector<vector<int>> result;
	cout << "start reading userInfoFile" << endl;
	for (; getline(userInfoFilePtr, temp);)
	{
		vector<int> rateTuple;
		vector<string> temp1 = split(temp, "\t");
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		int rate = std::atoi(temp1[2].c_str());
		rateTuple.push_back(uidx);
		rateTuple.push_back(iidx);
		rateTuple.push_back(rate);
		result.push_back(rateTuple);
	}
	userInfoFilePtr.close();
	return result;
}
map<int, map<int, int>> statisticRateitem(string filename)
{
	ifstream userInfoFilePtr;
	userInfoFilePtr.open(filename);
	string temp;
	map<int, map<int, int>> result;
	map<int, map<int, int>>::iterator itr;
	cout << "start reading userInfoFile" << endl;
	for (; getline(userInfoFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		int rate = std::atoi(temp1[2].c_str());
		itr = result.find(iidx);
		if (itr != result.end())
		{
			itr->second.insert(pair<int, int>(uidx, rate));
		}
		else
		{
			map<int, int> temp;
			temp.insert(pair<int, int>(uidx, rate));
			result.insert(pair<int, map<int, int>>(iidx, temp));
		}
	}
	userInfoFilePtr.close();
	return result;
}
double similarity(map<int, int> userArate, map<int, int> userBrate)
{
	map<int, int>::iterator itrA = userArate.begin();
	map<int, int>::iterator itrB = userBrate.begin();
	vector<int> intersection;
	double avgA = 0.0, avgB = 0.0,aa=0.0,bb=0.0,cc=0.0;
	while (itrA != userArate.end() && itrB != userBrate.end())
	{
		if (itrA->first == itrB->first)
		{
			intersection.push_back(itrA->first);
			itrA++;
			itrB++;
		}
		else if (itrA->first > itrB->first)
		{
			itrB++;
		}
		else
			itrA++;		
	}
	if (intersection.size() > 0)
	{
		for (int k = 0; k < intersection.size(); k++)
		{
			avgA += userArate[intersection[k]];
			avgB += userBrate[intersection[k]];
		}
		avgA = avgA / intersection.size();
		avgB = avgB / intersection.size();
		for (int k = 0; k < intersection.size(); k++)
		{
			aa += (userArate[intersection[k]] - avgA)*(userBrate[intersection[k]] - avgB);
			bb += (userArate[intersection[k]] - avgA)*(userArate[intersection[k]] - avgA);
			cc += (userBrate[intersection[k]] - avgB)*(userBrate[intersection[k]] - avgB);
		}
		if (bb == 0.0 || cc == 0.0)
		{
			aa = 0.0;
			bb = 0.0;
			cc = 0.0;
			for (int k = 0; k < intersection.size(); k++)
			{
				aa += userArate[intersection[k]] * userBrate[intersection[k]];
				bb += userArate[intersection[k]] * userArate[intersection[k]];
				cc += userBrate[intersection[k]] * userBrate[intersection[k]];
			}
		}
		bb = sqrt(bb);
		cc = sqrt(cc);
		return  aa / (bb*cc);
	}
	else
		return 0;
}
map<int, double> calculateAllSim(map<int, map<int, int>> rate, int user)
{
	map<int, double> result;
	map<int, map<int, int>>::iterator itr = rate.begin();
	for (; itr != rate.end(); itr++)
	{
		if (itr->first != user)
		{
			result.insert(pair<int, double>(itr->first, similarity(rate[user], itr->second)));
		}
	}
	return result;
}
int cmp(const PAIR &x, const PAIR &y)
{
	return x.second > y.second;
}
int cmp1(const PAIR &x, const PAIR &y)
{
	return x.second < y.second;
}
map<int, map<int,double>> findTopN(map<int, double> simS, int N)//map[0]表示最相似的N个，map[1]表示最不相似的N个
{
	map<int,double> top,down;
	vector<PAIR> pair_vec;
	for (map<int, double>::iterator map_iter = simS.begin(); map_iter != simS.end(); map_iter++)
	{
		pair_vec.push_back(make_pair(map_iter->first, map_iter->second));
	}
	sort(pair_vec.begin(), pair_vec.end(), cmp);
	int count = 0;
	for (vector<PAIR>::iterator curr = pair_vec.begin(); curr != pair_vec.end(); curr++)
	{
		
		
		if (count < N)
		{
			top.insert(pair<int, double>(curr->first, curr->second));
		}
		else if (count>=pair_vec.size()-N)
		{
			down.insert(pair<int, double>(curr->first, curr->second));
		}
		count++;
	}
	map<int, map<int,double>> result;
	result[0] = top;
	result[1] = down;
	return result;
}
map<int, map<int, double>> findTopN(map<int, double> simS, int N,double thread)//map[0]表示最相似的N个，map[1]表示最不相似的N个
{
	map<int, double> top, down;
	vector<PAIR> pair_vec;
	for (map<int, double>::iterator map_iter = simS.begin(); map_iter != simS.end(); map_iter++)
	{
		pair_vec.push_back(make_pair(map_iter->first, map_iter->second));
	}
	sort(pair_vec.begin(), pair_vec.end(), cmp);
	for (vector<PAIR>::iterator curr = pair_vec.begin(); curr != pair_vec.end() && top.size()< N; curr++)
	{
		if (curr->second>=thread)
		{
			top.insert(pair<int, double>(curr->first, curr->second));
		}
		else
			break;
	}
	sort(pair_vec.begin(), pair_vec.end(), cmp1);
	for (vector<PAIR>::iterator curr = pair_vec.begin(); curr != pair_vec.end() && down.size()< N; curr++)
	{
		if (curr->second<=(-1.0*thread))
		{
			down.insert(pair<int, double>(curr->first, curr->second));
		}
		else
			break;
	}
	map<int, map<int, double>> result;
	result[0] = top;
	result[1] = down;
	return result;
}

pair<map<int, map<int, double>>, map<int, map<int, double>>> findTopNall(map<int, map<int, int>> rateAll, int N)
{
	map<int, map<int, int>>::iterator itrRate = rateAll.begin();
	map<int, map<int, double>> topNall, downNall;
	for (; itrRate != rateAll.end(); itrRate++)
	{
		map<int, map<int, double>> result;
		int user = itrRate->first;
		result = findTopN(calculateAllSim(rateAll, user), N);
		topNall.insert(pair<int, map<int, double>>(user, result[0]));
		downNall.insert(pair<int, map<int, double>>(user, result[1]));
	}
	return pair<map<int, map<int, double>>, map<int, map<int, double>>>(topNall, downNall);
}
map<int, double> findTopN(user userInfo, int N, double thread)//map[0]表示最相似的N个，map[1]表示最不相似的N个
{
	int top = 0;
	map<int, double> result;
	vector<PAIR> pair_vec;
	map<int, double>::iterator itr = userInfo.sim.begin();
	for (; itr !=userInfo.sim.end(); itr++)
	{
		pair_vec.push_back(pair<int, double>(itr->first, itr->second));
		
	}
	sort(pair_vec.begin(), pair_vec.end(), cmp);
	for (vector<PAIR>::iterator curr = pair_vec.begin(); top< N && curr != pair_vec.end(); curr++)
	{
		if (curr->second >= thread)
		{
			result[curr->first] = curr->second;
			top++;
		}
		else
			break;
	}
	return result;
}
map<int,double> findTopN(item itemInfo, int N, double thread)//map[0]表示最相似的N个，map[1]表示最不相似的N个
{
	int top = 0;
	map<int, double> result;
	vector<PAIR> pair_vec;
	map<int, double>::iterator itr = itemInfo.sim.begin();
	for (; itr != itemInfo.sim.end(); itr++)
	{
		pair_vec.push_back(pair<int, double>(itr->first, itr->second));

	}
	sort(pair_vec.begin(), pair_vec.end(), cmp);
	for (vector<PAIR>::iterator curr = pair_vec.begin(); top< N && curr != pair_vec.end(); curr++)
	{
		if (curr->second >= thread)
		{
			result[curr->first] = curr->second;
			top++;
		}
		else
			break;
	}
	return result;
}
map<int, map<int, double>> readSimFile(string simFileName)
{
	map<int, map<int, double>> result;
	ifstream simFilePtr;
	simFilePtr.open(simFileName);
	string temp;
	for (; getline(simFilePtr, temp);)
	{
		map<int, double> temp2;
		vector<string> temp1 = split(temp, "\t");
		int idx = std::atoi(temp1[0].c_str());
		for (int i = 1; i < temp1.size(); i++)
		{
			vector<string> temp3 = split(temp1[i], ":");
			temp2.insert(pair<int, double>(std::atoi(temp3[0].c_str()), std::atof(temp3[1].c_str())));
		}
		result.insert(pair<int, map<int, double>>(idx, temp2));
	}
	simFilePtr.close();
	return result;
}
double vectorP(int a[], int b[], int n)
{
	double product = 0, sum1 = 0, sum2 = 0 ;
	for (int i = 0; i < n; i++)
	{
		product += a[i] * b[i];
		sum1 += a[i] * a[i];
		sum2 += b[i] * b[i];
	}
	if (sum1 == 0||sum2 == 0)
		product = 0;
	else
		product /= (sqrt(sum1)*sqrt(sum2));
	return product;
}

void findTopN(user userInfo, int N, double thread, double alpha ,double beta)//map[0]表示最相似的N个，map[1]表示最不相似的N个
{
	int top=0;
	vector<PAIR> pair_vec;
	for (int i = 1; i < userInfo.simA.size();i++)
	{
		if (i != userInfo.id)
		{
			double sim = alpha * userInfo.simS[i] + beta * userInfo.simA[i];
			userInfo.sim[i] = sim;
			pair_vec.push_back(pair<int,double> (i,sim));
		}
	}	
	sort(pair_vec.begin(), pair_vec.end(), cmp);
	for (vector<PAIR>::iterator curr = pair_vec.begin(); top< N && curr != pair_vec.end(); curr++)
	{
		if (curr->second >= thread)
		{
			userInfo.topN[curr->first] = curr->second;
			top++;
		}
		else
			break;
	}
}
void findTopN(item itemInfo, int N, double thread, double alpha, double beta)//map[0]表示最相似的N个，map[1]表示最不相似的N个
{
	int top=0;
	vector<PAIR> pair_vec;
	for (int i = 1; i < itemInfo.simA.size(); i++)
	{
		if (i != itemInfo.id)
		{
			double sim = alpha * itemInfo.simS[i] + beta * itemInfo.simA[i];
			itemInfo.sim[i] = sim;
			pair_vec.push_back(pair<int, double>(i, sim));
		}
	}
	sort(pair_vec.begin(), pair_vec.end(), cmp);
	for (vector<PAIR>::iterator curr = pair_vec.begin(); top< N && curr != pair_vec.end(); curr++)
	{
		if (curr->second >= thread)
		{
			itemInfo.topN[curr->first] = curr->second;
			top++;
		}
		else
			break;
	}
}