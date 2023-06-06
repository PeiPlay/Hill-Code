#pragma once
#include <vector>
#include <string>
#include <iostream>
#include "Matrix.h"
#define NUM_PER_UNIT 4

class CiphertextMatrix : public Matrix<double>//密文矩阵
{
public:
	CiphertextMatrix(const Matrix<double>& matrix, int npu) :Matrix<double>(matrix), num_per_unit(npu) {}
	CiphertextMatrix(uint32_t m, uint32_t n, int npu) : Matrix<double>(m, n), num_per_unit(npu) {}
	CiphertextMatrix(const CiphertextMatrix& matrix) : Matrix<double>(matrix), num_per_unit(matrix.num_per_unit) {}
	int& Spaces(){ return space; }
	int GetNumPerUnit() const { return num_per_unit; }
	std::string DecToString(const Matrix<double>& key)
	{
		Matrix<double> dec = key * (*this);
		int strlen = dec.GetRowNum() * num_per_unit;
		std::string result;
		for (int i = dec.GetRowNum() - 1; i >= 0; i--)
		{
			uint64_t bit_info = (dec.getUnit(i, 0) + 0.5);
			for (int j = 0; j < num_per_unit; j++)
			{
				result.push_back(bit_info % 1000);
				bit_info /= 1000;
			}
		}
		//字符串反转
		for (int i = 0; i < strlen / 2; i++)
		{
			char temp = result[i];
			result[i] = result[strlen - i - 1];
			result[strlen - i - 1] = temp;
		}

		//删除末尾多余的空格
		for (int i = 0; i < space; i++)
		{
			result.pop_back();
		}
		return result;
	}
private:
	int num_per_unit;//每个元素储存的数字个数
	int space;//添加的空格数
};

class OriInfo
{
	struct unit
	{
		int num_per_unit;
		std::string message;
	};
public:
	OriInfo(std::string str, int order)//str为输入的信息，order为矩阵的阶数
	{
		//对于一个order为n的密钥矩阵，与他相乘的信息矩阵为n*1的矩阵
		//这个n*1的矩阵每一个元素为double类型，它们每个最多储存3位数字
		//那么对于一个n阶密钥矩阵，它最多可以加密3n位的信息
		//所以将str分成多份，每一份为3n位，最后的一份不足3n位的需要补齐空格

		int length = str.size();
		int j = 0;
		for (int i = 0; i < length; i++)
		{
			if (i % (NUM_PER_UNIT * order) == 0)//每3n位分成一份
			{
				info.push_back(unit());
				j++;
				info[j - 1].num_per_unit = NUM_PER_UNIT;
			}
			info[j - 1].message.push_back(str[i]);
		}
		for (int i = length; i % (NUM_PER_UNIT * order) != 0; i++)
		{
			info[j - 1].message.push_back(' ');
		}
	}
	std::string GetInfo(int index)
	{
		return info[index].message;
	}
	int GetNumPerUnit(int index)
	{
		return info[index].num_per_unit;
	}
	int GetLength()
	{
		return info.size();
	}
	CiphertextMatrix ToMatrix_Enc(int index, const Matrix<double>& key) const
	{
		int length = info[index].message.size();
		int num_per_unit = info[index].num_per_unit;
		int order = key.GetColNum();
		Matrix<double> tmp(order, 1);
		
		int k = 0;
		int space = 0;
		for (int i = 0; i < order; i++)
		{
			for (int j = 0; j < num_per_unit; j++)
			{
				tmp(i, 0) *= 1000;
				if (k < length)
				{
					tmp(i, 0) += info[index].message[k];
				}
				else
				{
					tmp(i, 0) += ' ';
					space += 1;
				}
				k++;
			}
		}
		CiphertextMatrix result(key * tmp, info[index].num_per_unit);
		result.Spaces() = space;
		return result;
	}
	void ShowInfo() const
	{
		for (int i = 0; i < info.size(); i++)
		{
			std::cout << info[i].message<< "/" << info[i].num_per_unit << std::endl;
		}
		std::cout << "*" << std::endl;
	}
private:
	std::vector<unit> info;
};


