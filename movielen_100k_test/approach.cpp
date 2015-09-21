#include "similarity.h"
#include "classObject.h"
#include <time.h>

int RSVD(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, double lr=0.01, double hp1=0.1,int maxItr=10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i=0; i<train.size();i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor).transpose()* itemInfo[iidx].factor;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*userInfo[uidx].factor - hp1*itemInfo[iidx].factor);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor).transpose()* itemInfo[iidx].factor;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itemResult.txt"), resultFile("result.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
int BiasMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, double AVG,double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor).transpose()* itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp1*userInfo[uidx].factor);
			userInfo[uidx].bias += lr*(erro - hp2*userInfo[uidx].bias);
			itemInfo[iidx].factor += lr*(erro*userInfo[uidx].factor - hp1*itemInfo[iidx].factor);
			itemInfo[iidx].bias += lr*(erro - hp2*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor).transpose()* itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itemResult.txt"), resultFile("result.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		userResultFile << "\t" << itUser->second.bias;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		userResultFile << "\t" << itItem->second.bias;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
int SocialMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1,int factorNumber=20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor).transpose()* itemInfo[iidx].factor;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();
			
			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}
		
			itrSim = itemInfo[iidx].topN.begin();
			
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
		
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp2*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*userInfo[uidx].factor - hp2*itermRelationFactor - hp1*itemInfo[iidx].factor);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor).transpose()* itemInfo[iidx].factor;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itemResult.txt"), resultFile("result.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
int AttributeMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test,age age_case,sex sex_case,occupation occupation_case, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itermResult.txt"), resultFile("result.txt"), ageResultFile("ageResult.txt"), sexResultFile("sexResult.txt"), occupationResultFile("occupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
int PSMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1,int factorNumber=20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();
			
			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}
		
			itrSim = itemInfo[iidx].topN.begin();
			
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
			
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor;
		double erro = fabs( rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itermResult.txt"), resultFile("result.txt"), ageResultFile("ageResult.txt"), sexResultFile("sexResult.txt"), occupationResultFile("occupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
int BiasNASOMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();
			
			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}
		
			itrSim = itemInfo[iidx].topN.begin();
			
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
		
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itermResult.txt"), resultFile("result.txt"), ageResultFile("ageResult.txt"), sexResultFile("sexResult.txt"), occupationResultFile("occupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}

int BiasNSOMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();

			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}

			itrSim = itemInfo[iidx].topN.begin();

			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}

			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itermResult.txt"), resultFile("result.txt"), ageResultFile("ageResult.txt"), sexResultFile("sexResult.txt"), occupationResultFile("occupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
int BiasNAOMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();

			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}

			itrSim = itemInfo[iidx].topN.begin();

			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}

			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itermResult.txt"), resultFile("result.txt"), ageResultFile("ageResult.txt"), sexResultFile("sexResult.txt"), occupationResultFile("occupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}

int BiasNASMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();

			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}

			itrSim = itemInfo[iidx].topN.begin();

			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}

			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("userResult.txt"), itemResultFile("itermResult.txt"), resultFile("result.txt"), ageResultFile("ageResult.txt"), sexResultFile("sexResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}


int BiasNSMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();

			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}

			itrSim = itemInfo[iidx].topN.begin();

			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}

			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream resultFile("BiasNSMFresult.txt");
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
double GlobalMean(vector<vector<int>> train, vector<vector<int>>test)
{
	clock_t start, finish;
	start = clock();
	double sumRate = 0.0;
	for (int i = 0; i<train.size(); i++)
	{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			sumRate += rate;
	}
	double avgRate = sumRate / train.size();
	finish = clock();
	double dur = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "	total time : " << dur << std::endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double erro =fabs( rate - avgRate);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream resultFile("GlobalMeanresult.txt");
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "avgRate:" << avgRate << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return avgRate;
}
double UserMean(map<int, user> userInfo, vector<vector<int>> train, vector<vector<int>>test)
{
	clock_t start, finish;
	start = clock();
	map<int,int> sumRateCount;
	for (int i = 0; i<train.size(); i++)
	{
		int uidx = train[i][0];
		int iidx = train[i][1];
		int rate = train[i][2];
		userInfo[uidx].mean += rate;
		if (sumRateCount.find(uidx) != sumRateCount.end())
			sumRateCount[uidx] += 1;
		else
			sumRateCount.insert(pair<int, int>(uidx, 1));
	}
	for (int i = 1; i <= userInfo.size(); i++)
	{
		userInfo[i].mean = userInfo[i].mean / sumRateCount[i];
	}
	finish = clock();
	double dur = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "	total time : " << dur << std::endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double erro = fabs(rate - userInfo[uidx].mean);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream resultFile("UserMeanResult.txt"),userResultFile("UserMeanUserAvgResult.text");
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile.close();
	for (int i = 1; i <= userInfo.size(); i++)
	{
		userResultFile << i << "\t" << userInfo[i].mean << endl;
	}
	userResultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return MAP;
}

int BiasASOMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);

			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasASOMFuserResult.txt"), itemResultFile("BiasASOMFitermResult.txt"), resultFile("BiasASOMFresult.txt"), ageResultFile("BiasASOMFageResult.txt"), sexResultFile("BiasASOMFsexResult.txt"), occupationResultFile("BiasASOMFoccupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}
double ItemMean(map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test)
{
	clock_t start, finish;
	start = clock();
	map<int, int> sumRateCount;
	for (int i = 0; i<train.size(); i++)
	{
		int uidx = train[i][0];
		int iidx = train[i][1];
		int rate = train[i][2];
		itemInfo[iidx].mean += rate;
		if (sumRateCount.find(iidx) != sumRateCount.end())
			sumRateCount[iidx] += 1;
		else
			sumRateCount.insert(pair<int, int>(iidx, 1));
	}
	for (int i = 1; i <= itemInfo.size(); i++)
	{
		if (sumRateCount.find(i) != sumRateCount.end())
			itemInfo[i].mean /= sumRateCount[i];
	}
	finish = clock();
	double dur = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "	total time : " << dur << std::endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		if (itemInfo[iidx].mean != 0)
		{
			double erro = fabs(rate - itemInfo[iidx].mean);
			testCount++;
			MAP = (MAP*(testCount - 1) + erro) / testCount;
			RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
		}			
	}
	ofstream resultFile("ItemMeanResult.txt"), itemResultFile("ItemMeanItemAvgResult.text");
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile.close();
	for (int i = 1; i <= itemInfo.size(); i++)
	{
		itemResultFile << i << "\t" << itemInfo[i].mean << endl;
	}
	itemResultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return MAP;
}

int BiasNAMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();

			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}

			itrSim = itemInfo[iidx].topN.begin();
	
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
		
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] ) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasNAMFUserResult.txt"), itemResultFile("BiasNAMFItermResult.txt"), resultFile("BiasNAMFResult.txt"), ageResultFile("BiasNAMFAgeResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}


int BiasNOMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	int itrTime = 0;
	clock_t start, finish;
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();
			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}

			itrSim = itemInfo[iidx].topN.begin();
			
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
		
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasNOMFuserResult.txt"), itemResultFile("BiasNOMFitermResult.txt"), resultFile("BiasNOMFresult.txt"),occupationResultFile("BiasNOMFoccupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}

int BiasNASOPMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	double beta = 0.1;
	int itrTime = 0;
	clock_t start, finish;
	double max = 0;
	for (int i = 0; i < train.size(); i++)
	{
		int uidx = train[i][0];
		int iidx = train[i][1];
		int rate = train[i][2];
		itemInfo[iidx].popular += 1;
		if (itemInfo[iidx].popular>max)
			max = itemInfo[iidx].popular;
	}
	for (int i = 1; i <= itemInfo.size(); i++)
	{
		itemInfo[i].popular /= max;
	}
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular ;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();
	
			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}
			
			itrSim = itemInfo[iidx].topN.begin();
			
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
		
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
			beta += lr*(erro *itemInfo[iidx].popular- hp4*beta);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasNASOPMFuserResult.txt"), itemResultFile("BiasNASOPMFitermResult.txt"), resultFile("BiasNASOPMFresult.txt"), ageResultFile("BiasNASOPMFageResult.txt"), sexResultFile("BiasNASOPMFsexResult.txt"), occupationResultFile("BiasNASOPMFoccupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile << "Popular Weight:" << beta << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}




int BiasNASOPCMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	double beta = 0.1;
	int itrTime = 0;
	clock_t start, finish;
	double max = 0;
	catigries catigries_case;
	double catitroyCount[19] = { 0 };
	for (int i = 0; i < train.size(); i++)
	{
		int uidx = train[i][0];
		int iidx = train[i][1];
		int rate = train[i][2];
		itemInfo[iidx].popular += 1;
		for (int k = 0; k<19; k++)
		{
			if (itemInfo[iidx].catigories[k] == 1)
			{
				catitroyCount[k] += 1;
			}		
		}
		if (itemInfo[iidx].popular>max)
			max = itemInfo[iidx].popular;	
	}
	double maxCatigries = catitroyCount[0];
	for (int k = 1; k<19; k++)
	{
		if (catitroyCount[k] > maxCatigries)
		{
			maxCatigries=catitroyCount[k];
		}
	}
	for (int k = 0; k<19; k++)
	{
		catigries_case.popular[k] = catitroyCount[k] / maxCatigries;
	}
	for (int i = 1; i <= itemInfo.size(); i++)
	{
		itemInfo[i].popular /= max;
	}
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			vector<int> cat;
			double catPopular = 0.0;
			for (int k = 0; k < 19; k++)
			{
				if (itemInfo[iidx].catigories[k] == 1)
				{
					cat.push_back(k);
					catPopular += catigries_case.popular[k] * catigries_case.factor[k];
				}			
			}
			if (cat.size()>0)
				catPopular /= cat.size();
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular + catPopular;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
			VectorXd socialFactor = VectorXd::Zero(factorNumber), itermRelationFactor = VectorXd::Zero(factorNumber);
			map<int, double>::iterator itrSim = userInfo[uidx].topN.begin();
		
			for (; itrSim != userInfo[uidx].topN.end(); itrSim++)
			{
				socialFactor += itrSim->second*(userInfo[uidx].factor - userInfo[itrSim->first].factor);
			}
		
			itrSim = itemInfo[iidx].topN.begin();
		
			for (; itrSim != itemInfo[iidx].topN.end(); itrSim++)
			{
				itermRelationFactor += itrSim->second*(itemInfo[iidx].factor - itemInfo[itrSim->first].factor);
			}
	
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor - hp3*socialFactor - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]) - hp3* itermRelationFactor - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
			beta += lr*(erro *itemInfo[iidx].popular - hp4*beta);
			for (int k = 0; k < cat.size(); k++)
			{
				catigries_case.factor[cat[k]] += lr*(erro *catigries_case.popular[k] - hp4*catigries_case.factor[cat[k]]);
			}
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		vector<int> cat;
		double catPopular = 0.0;
		for (int k = 0; k < 19; k++)
		{
			if (itemInfo[iidx].catigories[k] == 1)
			{
				cat.push_back(k);
				catPopular += catigries_case.popular[k] * catigries_case.factor[k];
			}
		}
		if (cat.size()>0)
			catPopular /= cat.size();
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular + catPopular;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasNASOPCMFuserResult.txt"), itemResultFile("BiasNASOPCMFitermResult.txt"), resultFile("BiasNASOPCMFresult.txt"), ageResultFile("BiasNASOPCMFageResult.txt"), sexResultFile("BiasNASOPCMFsexResult.txt"), occupationResultFile("BiasNASOPCMFoccupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile << "Popular Weight:" << beta << endl;
	for (int k = 0; k < 19; k++)
	{
		resultFile<<"Catigory"<< k <<"Popular Weight:"<<catigries_case.factor[k]<<endl;
	}
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}

int BiasASOPMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	double beta = 0.1;
	int itrTime = 0;
	clock_t start, finish;
	double max = 0;
	for (int i = 0; i < train.size(); i++)
	{
		int uidx = train[i][0];
		int iidx = train[i][1];
		int rate = train[i][2];
		itemInfo[iidx].popular += 1;
		if (itemInfo[iidx].popular>max)
			max = itemInfo[iidx].popular;
	}
	for (int i = 1; i <= itemInfo.size(); i++)
	{
		itemInfo[i].popular /= max;
	}
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
				
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor  - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]])  - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
			beta += lr*(erro *itemInfo[iidx].popular - hp4*beta);
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasASOPMFuserResult.txt"), itemResultFile("BiasASOPMFitermResult.txt"), resultFile("BiasASOPMFresult.txt"), ageResultFile("BiasASOPMFageResult.txt"), sexResultFile("BiasASOPMFsexResult.txt"), occupationResultFile("BiasASOPMFoccupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile << "Popular Weight:" << beta << endl;
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;
}


int BiasASOPCMF(map<int, user> userInfo, map<int, item> itemInfo, vector<vector<int>> train, vector<vector<int>>test, age age_case, sex sex_case, occupation occupation_case, double AVG, double lr = 0.01, double hp1 = 0.1, double hp2 = 0.1, double hp3 = 0.1, double hp4 = 0.1, int factorNumber = 20, int maxItr = 10000)
{
	double beta = 0.1;
	int itrTime = 0;
	clock_t start, finish;
	double max = 0;
	catigries catigries_case;
	double catitroyCount[19] = { 0 };
	for (int i = 0; i < train.size(); i++)
	{
		int uidx = train[i][0];
		int iidx = train[i][1];
		int rate = train[i][2];
		itemInfo[iidx].popular += 1;
		for (int k = 0; k<19; k++)
		{
			if (itemInfo[iidx].catigories[k] == 1)
			{
				catitroyCount[k] += 1;
			}
		}
		if (itemInfo[iidx].popular>max)
			max = itemInfo[iidx].popular;
	}
	double maxCatigries = catitroyCount[0];
	for (int k = 1; k<19; k++)
	{
		if (catitroyCount[k] > maxCatigries)
		{
			maxCatigries = catitroyCount[k];
		}
	}
	for (int k = 0; k<19; k++)
	{
		catigries_case.popular[k] = catitroyCount[k] / maxCatigries;
	}
	for (int i = 1; i <= itemInfo.size(); i++)
	{
		itemInfo[i].popular /= max;
	}
	double lastMAE = 1.0, thisMAE = 0.0;
	for (; itrTime < maxItr&&fabs(lastMAE - thisMAE)>0.000001; itrTime++)
	{
		start = clock();
		lastMAE = thisMAE;
		thisMAE = 0.0;
		for (int i = 0; i<train.size(); i++)
		{
			int uidx = train[i][0];
			int iidx = train[i][1];
			int rate = train[i][2];
			vector<int> cat;
			double catPopular = 0.0;
			for (int k = 0; k < 19; k++)
			{
				if (itemInfo[iidx].catigories[k] == 1)
				{
					cat.push_back(k);
					catPopular += catigries_case.popular[k] * catigries_case.factor[k];
				}
			}
			if (cat.size()>0)
				catPopular /= cat.size();
			double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular + catPopular;
			double erro = rate - ratePre;
			thisMAE += fabs(erro);
		
			userInfo[uidx].factor += lr*(erro*itemInfo[iidx].factor  - hp1*userInfo[uidx].factor);
			itemInfo[iidx].factor += lr*(erro*(userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]])  - hp1*itemInfo[iidx].factor);
			age_case.factor[userInfo[uidx].age / 25] += lr*(erro*itemInfo[iidx].factor - hp2*age_case.factor[userInfo[uidx].age / 25]);
			sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] += lr*(erro*itemInfo[iidx].factor - hp2*sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]]);
			occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]] += lr*(erro*itemInfo[iidx].factor - hp2*occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]);
			userInfo[uidx].bias += lr*(erro - hp4*userInfo[uidx].bias);
			itemInfo[iidx].bias += lr*(erro - hp4*itemInfo[iidx].bias);
			beta += lr*(erro *itemInfo[iidx].popular - hp4*beta);
			for (int k = 0; k < cat.size(); k++)
			{
				catigries_case.factor[cat[k]] += lr*(erro *catigries_case.popular[k] - hp4*catigries_case.factor[cat[k]]);
			}
		}
		if (itrTime % 10 == 0)
			lr *= 0.95;
		thisMAE /= train.size();
		finish = clock();
		double dur = (double)(finish - start) / CLOCKS_PER_SEC;
		std::cout << "loop:" << itrTime << "	total time : " << dur << "	thisMAE:" << thisMAE << std::endl;
	}
	cout << "iterator  over,start test" << endl;
	double RSME = 0.0, MAP = 0.0;
	int testCount = 0;
	for (int i = 0; i<test.size(); i++)
	{
		int uidx = test[i][0];
		int iidx = test[i][1];
		int rate = test[i][2];
		vector<int> cat;
		double catPopular = 0.0;
		for (int k = 0; k < 19; k++)
		{
			if (itemInfo[iidx].catigories[k] == 1)
			{
				cat.push_back(k);
				catPopular += catigries_case.popular[k] * catigries_case.factor[k];
			}
		}
		if (cat.size()>0)
			catPopular /= cat.size();
		double ratePre = (userInfo[uidx].factor + age_case.factor[userInfo[uidx].age / 25] + sex_case.factor[sex_case.nameToId[userInfo[uidx].sex]] + occupation_case.factor[occupation_case.nameToId[userInfo[uidx].occupation]]).transpose() * itemInfo[iidx].factor + AVG + userInfo[uidx].bias + itemInfo[iidx].bias + beta*itemInfo[iidx].popular + catPopular;
		double erro = fabs(rate - ratePre);
		testCount++;
		MAP = (MAP*(testCount - 1) + erro) / testCount;
		RSME = (RSME*(testCount - 1) + erro*erro) / testCount;
	}
	ofstream userResultFile("BiasASOPCMFuserResult.txt"), itemResultFile("BiasASOPCMFitermResult.txt"), resultFile("BiasASOPCMFresult.txt"), ageResultFile("BiasASOPCMFageResult.txt"), sexResultFile("BiasASOPCMFsexResult.txt"), occupationResultFile("BiasASOPCMFoccupationResult.txt");
	map<int, user>::iterator itUser;
	for (itUser = userInfo.begin(); itUser != userInfo.end(); itUser++)
	{
		userResultFile << itUser->first;
		for (int k = 0; k < itUser->second.factor.size(); k++)
		{
			userResultFile << "\t" << itUser->second.factor[k];
		}
		userResultFile << endl;
	}
	userResultFile.close();
	map<int, item>::iterator itItem;
	for (itItem = itemInfo.begin(); itItem != itemInfo.end(); itItem++)
	{
		itemResultFile << itItem->first;
		for (int k = 0; k < itItem->second.factor.size(); k++)
		{
			itemResultFile << "\t" << itItem->second.factor[k];
		}
		itemResultFile << endl;
	}
	itemResultFile.close();
	map<int, VectorXd>::iterator itAge;
	for (itAge = age_case.factor.begin(); itAge != age_case.factor.end(); itAge++)
	{
		ageResultFile << itAge->first;
		for (int k = 0; k < itAge->second.size(); k++)
		{
			ageResultFile << "\t" << itAge->second[k];
		}
		ageResultFile << endl;
	}
	ageResultFile.close();
	map<int, VectorXd>::iterator itSex;
	for (itSex = sex_case.factor.begin(); itSex != sex_case.factor.end(); itSex++)
	{
		sexResultFile << itSex->first;
		for (int k = 0; k < itSex->second.size(); k++)
		{
			sexResultFile << "\t" << itSex->second[k];
		}
		sexResultFile << endl;
	}
	sexResultFile.close();
	map<int, VectorXd>::iterator itOccupation;
	for (itOccupation = occupation_case.factor.begin(); itOccupation != occupation_case.factor.end(); itOccupation++)
	{
		occupationResultFile << itOccupation->first;
		for (int k = 0; k < itOccupation->second.size(); k++)
		{
			occupationResultFile << "\t" << itOccupation->second[k];
		}
		occupationResultFile << endl;
	}
	occupationResultFile.close();
	resultFile << "MAP:" << MAP << endl;
	resultFile << "RSME:" << sqrt(RSME) << endl;
	resultFile << "Iterator time:" << itrTime << endl;
	resultFile << "Popular Weight:" << beta << endl;
	for (int k = 0; k < 19; k++)
	{
		resultFile << "Catigory" << k << "Popular Weight:" << catigries_case.factor[k] << endl;
	}
	resultFile.close();
	cout << MAP << endl;
	cout << sqrt(RSME);
	return itrTime;

}