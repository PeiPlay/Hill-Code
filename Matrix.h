#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <initializer_list>



//矩阵类，用于存储矩阵数据和进行矩阵运算
//矩阵为模板类，可以存储任意类型的数据

template <typename T>
class Matrix
{
public:
	Matrix(uint32_t m, uint32_t n);									//构造函数，m为行数，n为列数
	Matrix(const Matrix<T>& mat);									//拷贝构造函数
	Matrix(const T* mat, uint32_t m, uint32_t n);					//构造函数，mat为数据指针，m为行数，n为列数
	Matrix(std::initializer_list<std::initializer_list<T>> list);	//大括号初始化矩阵

	//运算符重载
	Matrix<T>& operator=(const Matrix<T>& mat);			//赋值运算符重载
	Matrix<T> operator+(const Matrix<T>& mat) const;	//加法运算符重载
	Matrix<T> operator-(const Matrix<T>& mat) const;	//减法运算符重载
	Matrix<T> operator*(const Matrix<T>& mat) const;	//乘法运算符重载
	Matrix<T> operator*(const T& num) const;			//数乘运算符重载
	Matrix<T>& operator+=(const Matrix<T>& mat);		//加法赋值运算符重载
	Matrix<T>& operator-=(const Matrix<T>& mat);		//减法赋值运算符重载
	Matrix<T>& operator*=(const Matrix<T>& mat);		//乘法赋值运算符重载
	Matrix<T>& operator*=(const T& num);				//数乘赋值运算符重载
	bool operator==(const Matrix<T>& mat) const;		//相等运算符重载
	bool operator!=(const Matrix<T>& mat) const;		//不等运算符重载
	

	//矩阵相关运算
	void Add(const Matrix<T>& mat);		//矩阵加法
	void Sub(const Matrix<T>& mat);		//矩阵减法
	void Mul(const Matrix<T>& mat);		//矩阵乘法
	void Mul(const T& num);				//矩阵数乘

	void setUnit(uint32_t m, uint32_t n, T u)
	{
		data[m][n] = u;
	}

	T getUnit(uint32_t m, uint32_t n) const
	{
		return data[m][n];
	}

	friend void GetRandomInvertibleMatrix(uint32_t seed)	//生成随机可逆矩阵
	{

	}

	void RowExchange(uint32_t i, uint32_t j);				//矩阵行交换
	void ColExchange(uint32_t i, uint32_t j);				//矩阵列交换
	void RowAdd(uint32_t i, uint32_t j, T k);				//矩阵行倍加
	void ColAdd(uint32_t i, uint32_t j, T k);				//矩阵列倍加
	

	//矩阵数据访问
	T& operator()(uint32_t i, uint32_t j);					//访问矩阵第i行第j列的数据
	const T& operator()(uint32_t i, uint32_t j) const;		//访问矩阵第i行第j列的数据
	Matrix<T> GetRow(uint32_t i) const;						//获取矩阵第i行
	Matrix<T> GetCol(uint32_t j) const;						//获取矩阵第j列
	uint32_t getRank() const;								//获取矩阵秩
	uint32_t getNullity() const;							//获取矩阵零度
	T getTrace() const;										//获取矩阵迹
	T getNorm() const;										//获取矩阵范数


	//矩阵相关操作
	Matrix<T> GetEchelonForm() const;			//获取矩阵最简形式
	Matrix<T>& ToEchelonForm();					//矩阵化为最简形式
	Matrix<T> GetReducedEchelonForm() const;	//获取矩阵行最简形式
	Matrix<T>& ToReducedEchelonForm();			//矩阵化为行最简形式


	Matrix<T>& Transpose();						//矩阵转置
	Matrix<T> GetTranspose() const;				//获取矩阵转置

	Matrix<T>& Inverse();						//矩阵求逆
	Matrix<T> GetInverse() const;				//获取矩阵逆矩阵
	bool isInvertible() const;					//判断矩阵是否可逆

	T Cofactor(uint32_t i, uint32_t j) const;	//矩阵求余子式

	Matrix<T>& Cofactorize();					//矩阵求余子式矩阵
	Matrix<T> GetCofactorization() const;		//获取矩阵余子式矩阵

	Matrix<T>& Adjoint();						//矩阵求伴随矩阵
	Matrix<T> GetAdjoint() const;				//获取矩阵伴随矩阵

	T Determinant() const;						//矩阵求行列式

	//矩阵相关信息
	void Print() const;							//打印矩阵
	void ShowInfo() const;						//显示矩阵信息

	uint32_t GetRowNum() const;					//获取矩阵行数
	uint32_t GetColNum() const;					//获取矩阵列数


private:
	T Abs(T num)
	{
		if (num >= 0)
		{
			return num;
		}
		else
		{
			return num * -1;
		}
	}

	//矩阵数据
	std::vector<std::vector<T>> data;
};

template <typename T>
Matrix<T>::Matrix(uint32_t m, uint32_t n)
{
	data.resize(m);
	for (uint32_t i = 0; i < m; i++)
	{
		data[i].resize(n);
	}
}
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& mat)
{
	data = mat.data;
}
template <typename T>
Matrix<T>::Matrix(const T* mat, uint32_t m, uint32_t n)
{
	data.resize(m);
	for (uint32_t i = 0; i < m; i++)
	{
		data[i].resize(n);
		for (uint32_t j = 0; j < n; j++)
		{
			data[i][j] = mat[i * n + j];
		}
	}
}
template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> list)
{
	//内侧的阵列代表行

	uint32_t m = list.size();
	uint32_t n = list.begin()->size();
	data.resize(m);
	uint32_t i = 0;
	for (auto& row : list)
	{
		data[i].resize(n);
		uint32_t j = 0;
		for (auto& col : row)
		{
			data[i][j] = col;
			j++;
		}
		i++;
	}
}
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat) const
{
	Matrix<T> res(*this);
	res.Add(mat);
	return res;
}
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& mat) const
{
	Matrix<T> res(*this);
	res.Sub(mat);
	return res;
}
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const
{
	Matrix<T> res(*this);
	res.Mul(mat);
	return res;
}
template <typename T>
Matrix<T> Matrix<T>::operator*(const T& num) const
{
	Matrix<T> res(*this);
	res.Mul(num);
	return res;
}
template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& mat)
{
	Add(mat);
	return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& mat)
{
	Sub(mat);
	return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& mat)
{
	Mul(mat);
	return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T& num)
{
	Mul(num);
	return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat)
{
	data = mat.data;
	return *this;
}
template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& mat) const
{
	if (data.size() != mat.data.size())
	{
		return false;
	}
	for (uint32_t i = 0; i < data.size(); i++)
	{
		if (data[i].size() != mat.data[i].size())
		{
			return false;
		}
		for (uint32_t j = 0; j < data[i].size(); j++)
		{
			if (data[i][j] != mat.data[i][j])
			{
				return false;
			}
		}
	}
	return true;
}
template <typename T>
bool Matrix<T>::operator!=(const Matrix<T>& mat) const
{
	return !(*this == mat);
}

template <typename T>
void Matrix<T>::RowExchange(uint32_t i, uint32_t j)				//矩阵行交换
{
	//交换行i和行j
	std::vector<T> temp = data[i];
	data[i] = data[j];
	data[j] = temp;
}
template <typename T>
void Matrix<T>::ColExchange(uint32_t i, uint32_t j)				//矩阵列交换
{
	//交换列i和列j
	for (uint32_t k = 0; k < data.size(); k++)
	{
		T temp = data[k][i];
		data[k][i] = data[k][j];
		data[k][j] = temp;
	}
}
template <typename T>
void Matrix<T>::RowAdd(uint32_t i, uint32_t j, T k)				//矩阵行倍加
{
	//将行j的k倍加到行i上
	for (uint32_t l = 0; l < data[i].size(); l++)
	{
		data[i][l] += data[j][l] * k;
	}
}
template <typename T>
void Matrix<T>::ColAdd(uint32_t i, uint32_t j, T k)				//矩阵列倍加
{
	//将列j的k倍加到列i上
	for (uint32_t l = 0; l < data.size(); l++)
	{
		data[l][i] += data[l][j] * k;
	}	
}


template <typename T>
uint32_t Matrix<T>::getRank() const
{
	Matrix<T> mat = GetEchelonForm();
	uint32_t m = mat.GetRowNum();
	uint32_t n = mat.GetColNum();
	uint32_t rank = 0;
	for (uint32_t i = 0; i < m; i++)
	{
		bool isZero = true;
		for (uint32_t j = 0; j < n; j++)
		{
			if (mat(i, j) != (T)0)
			{
				isZero = false;
				break;
			}
		}
		if (!isZero)
		{
			rank++;
		}
	}
	return rank;
}
template <typename T>
uint32_t Matrix<T>::getNullity() const
{
	return GetColNum() - getRank();
}
template <typename T>
T Matrix<T>::getTrace() const
{
	uint32_t m = GetRowNum();
	uint32_t n = GetColNum();
	T trace = (T)0;
	for (uint32_t i = 0; i < m && i < n; i++)
	{
		trace += data[i][i];
	}
	return trace;
}
template <typename T>
T Matrix<T>::getNorm() const
{
	uint32_t m = GetRowNum();
	uint32_t n = GetColNum();
	T norm = (T)0;
	for (uint32_t i = 0; i < m; i++)
	{
		for (uint32_t j = 0; j < n; j++)
		{
			norm += data[i][j] * data[i][j];
		}
	}
	return std::sqrt(norm);
}
template <typename T>
Matrix<T>& Matrix<T>::ToEchelonForm()
{
	uint32_t m = GetRowNum();
	uint32_t n = GetColNum();
	uint32_t i = 0;
	uint32_t j = 0;
	while (i < m && j < n)
	{
		//找到列j中第i行及其之后的最大元素所在的行
		uint32_t maxRow = i;
		for (uint32_t k = i + 1; k < m; k++)
		{

			if (Abs(data[k][j]) > Abs(data[maxRow][j]))
			{
				maxRow = k;
			}
		}
		//如果最大元素为0，则列j全为0，处理下一列
		if (data[maxRow][j] == 0)
		{
			j++;
		}
		else
		{
			//交换第i行和第maxRow行
			if (i != maxRow)
			{
				for (uint32_t k = j; k < n; k++)
				{
					std::swap(data[i][k], data[maxRow][k]);
				}
			}
			//将第i行的首元素缩放为1
			T scale = data[i][j];
			for (uint32_t k = j; k < n; k++)
			{
				data[i][k] /= scale;
			}
			//将第i行以下的第j列元素消为0
			for (uint32_t k = i + 1; k < m; k++)
			{
				T scale = data[k][j];
				for (uint32_t l = j; l < n; l++)
				{
					data[k][l] -= scale * data[i][l];
				}
			}
			i++;
			j++;
		}
	}
	return *this;
}
template <typename T>
Matrix<T> Matrix<T>::GetEchelonForm() const
{
	Matrix<T> res(*this);
	return res.ToEchelonForm();
}
template <typename T>
Matrix<T>& Matrix<T>::ToReducedEchelonForm()
{
	uint32_t m = GetRowNum();
	uint32_t n = GetColNum();
	uint32_t i = 0;
	uint32_t j = 0;
	while (i < m && j < n)
	{
		//找到列j中第i行及其之后的最大元素所在的行
		uint32_t maxRow = i;
		for (uint32_t k = i + 1; k < m; k++)
		{
			if (Abs(data[k][j]) > Abs(data[maxRow][j]))
			{
				maxRow = k;
			}
		}
		//如果最大元素为0，则列j全为0，处理下一列
		if (data[maxRow][j] == 0)
		{
			j++;
		}
		else
		{
			//交换第i行和第maxRow行
			if (i != maxRow)
			{
				for (uint32_t k = j; k < n; k++)
				{
					std::swap(data[i][k], data[maxRow][k]);
				}
			}
			//将第i行的首元素缩放为1
			T scale = data[i][j];
			for (uint32_t k = j; k < n; k++)
			{
				data[i][k] /= scale;
			}
			//将第i行以下的第j列元素消为0
			for (uint32_t k = i + 1; k < m; k++)
			{
				T scale = data[k][j];
				for (uint32_t l = j; l < n; l++)
				{
					data[k][l] -= scale * data[i][l];
				}
			}
			//将第i行以上的第j列元素消为0
			for (uint32_t k = 0; k < i; k++)
			{
				T scale = data[k][j];
				for (uint32_t l = j; l < n; l++)
				{
					data[k][l] -= scale * data[i][l];
				}
			}
			i++;
			j++;
		}
	}

	return *this;
}
template <typename T>
Matrix<T> Matrix<T>::GetReducedEchelonForm() const
{
	Matrix<T> res(*this);
	return res.ToReducedEchelonForm();
}
template <typename T>
uint32_t Matrix<T>::GetRowNum() const
{
	return data.size();
}
template <typename T>
uint32_t Matrix<T>::GetColNum() const
{
	return data[0].size();
}
template <typename T>
T& Matrix<T>::operator()(uint32_t i, uint32_t j)
{
	return data[i][j];
}
template <typename T>
const T& Matrix<T>::operator()(uint32_t i, uint32_t j) const
{
	return data[i][j];
}
template <typename T>
Matrix<T> Matrix<T>::GetRow(uint32_t i) const
{
	Matrix<T> res(1, GetColNum());
	for (uint32_t j = 0; j < GetColNum(); j++)
	{
		res(0, j) = data[i][j];
	}
	return res;
}
template <typename T>
Matrix<T> Matrix<T>::GetCol(uint32_t j) const
{
	Matrix<T> res(GetRowNum(), 1);
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		res(i, 0) = data[i][j];
	}
	return res;
}
template <typename T>
void Matrix<T>::Add(const Matrix<T>& mat)
{
	if (GetRowNum() != mat.GetRowNum() || GetColNum() != mat.GetColNum())
	{
		std::cout << "Error: Matrix::Add: Matrix size not match!" << std::endl;
		return;
	}
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			data[i][j] += mat(i, j);
		}
	}
}
template <typename T>
void Matrix<T>::Sub(const Matrix<T>& mat)
{
	if (GetRowNum() != mat.GetRowNum() || GetColNum() != mat.GetColNum())
	{
		std::cout << "Error: Matrix::Sub: Matrix size not match!" << std::endl;
		return;
	}
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			data[i][j] -= mat(i, j);
		}
	}
}
template <typename T>
void Matrix<T>::Mul(const Matrix<T>& mat)
{
	if (GetColNum() != mat.GetRowNum())
	{
		std::cout << "Error: Matrix::Mul: Matrix size not match!" << std::endl;
		return;
	}
	Matrix<T> res(GetRowNum(), mat.GetColNum());
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < mat.GetColNum(); j++)
		{
			for (uint32_t k = 0; k < GetColNum(); k++)
			{
				res(i, j) += data[i][k] * mat(k, j);
			}
		}
	}
	data = res.data;
}
template <typename T>
void Matrix<T>::Mul(const T& val)
{
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			data[i][j] *= val;
		}
	}
}
template <typename T>
Matrix<T>& Matrix<T>::Transpose()
{
	//实现矩阵转置
	Matrix<T> res(GetColNum(), GetRowNum());
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			res(j, i) = data[i][j];
		}
	}
	data = res.data;
	return *this;
}
template <typename T>
Matrix<T> Matrix<T>::GetTranspose() const
{
	Matrix<T> res(*this);
	return res.Transpose();
}
template <typename T>
T Matrix<T>::Determinant() const
{
	//实现矩阵求行列式
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::Determinant: Matrix size not match!" << std::endl;
		return false;
	}
	if (GetRowNum() == 1)
	{
		return data[0][0];
	}
	else if (GetRowNum() == 2)
	{
		return data[0][0] * data[1][1] - data[0][1] * data[1][0];
	}
	else
	{
		T res = 0;
		for (uint32_t i = 0; i < GetRowNum(); i++)
		{
			Matrix<T> tmp(GetRowNum() - 1, GetColNum() - 1);
			for (uint32_t j = 1; j < GetRowNum(); j++)
			{
				for (uint32_t k = 0; k < GetColNum(); k++)
				{
					if (k < i)
					{
						tmp(j - 1, k) = data[j][k];
					}
					else if (k > i)
					{
						tmp(j - 1, k - 1) = data[j][k];
					}
				}
			}
			res += data[0][i] * tmp.Determinant() * ((i % 2 == 0) ? 1 : -1);
		}
		return res;
	}
}
template <typename T>
bool Matrix<T>::isInvertible() const
{
	//实现判断矩阵是否可逆
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::isInvertible: Only square matrix can be inverted!" << std::endl;
		return false;
	}
	return Determinant() != 0;
}
template <typename T>
Matrix<T>& Matrix<T>::Inverse()
{
	//实现矩阵求逆
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::Inverse: Only square matrix can be inverted!" << std::endl;
		return *this;
	}
	//使用高斯消元法求逆
	Matrix<T> augmentation(GetRowNum(), GetColNum() * 2);	//增广矩阵
	//初始化增广矩阵
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			augmentation(i, j) = data[i][j];
		}
		augmentation(i, GetColNum() + i) = 1;
	}
	//高斯消元
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		//找到第i列的主元
		uint32_t max_row = i;
		for (uint32_t j = i + 1; j < GetRowNum(); j++)
		{
			if (augmentation(j, i) > augmentation(max_row, i))
			{
				max_row = j;
			}
		}
		//交换第i行和第max_row行
		for (uint32_t j = 0; j < GetColNum() * 2; j++)
		{
			std::swap(augmentation(i, j), augmentation(max_row, j));
		}
		//将第i行的主元缩放为1
		T scale = augmentation(i, i);
		for (uint32_t j = 0; j < GetColNum() * 2; j++)
		{
			augmentation(i, j) /= scale;
		}
		//将第i列的其他元素消为0
		for (uint32_t j = 0; j < GetRowNum(); j++)
		{
			if (j != i)
			{
				T scale = augmentation(j, i);
				for (uint32_t k = 0; k < GetColNum() * 2; k++)
				{
					augmentation(j, k) -= scale * augmentation(i, k);
				}
			}
		}
	}
	//将增广矩阵的右半部分作为结果
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			data[i][j] = augmentation(i, GetColNum() + j);
		}
	}
	return *this;
}
template <typename T>
Matrix<T> Matrix<T>::GetInverse() const
{
	//实现矩阵求逆
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::GetInverse: Only square matrix can be inverted!" << std::endl;
		return *this;
	}
	if (Determinant() == 0)
	{
		std::cout << "Error: Matrix::GetInverse: Matrix is not invertible!" << std::endl;
		return *this;
	}
	Matrix<T> res(GetRowNum(), GetColNum());
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			Matrix<T> tmp(GetRowNum() - 1, GetColNum() - 1);
			for (uint32_t k = 0; k < GetRowNum(); k++)
			{
				for (uint32_t l = 0; l < GetColNum(); l++)
				{
					if (k < i && l < j)
					{
						tmp(k, l) = data[k][l];
					}
					else if (k < i && l > j)
					{
						tmp(k, l - 1) = data[k][l];
					}
					else if (k > i && l < j)
					{
						tmp(k - 1, l) = data[k][l];
					}
					else if (k > i && l > j)
					{
						tmp(k - 1, l - 1) = data[k][l];
					}
				}
			}
			res(i, j) = tmp.Determinant() * ((i + j) % 2 == 0 ? (T)1 : (T)-1);
		}
	}
	res.Transpose();
	res.Mul((T)1 / Determinant());
	return res;
}
template <typename T>
void Matrix<T>::Print() const
{
	//打印矩阵
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			std::cout << data[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
template <typename T>
void Matrix<T>::ShowInfo() const
{
	//打印矩阵信息
	std::cout << "***************************************" << std::endl;
	std::cout << "Matrix Info:" << std::endl;
	std::cout << "RowNum: " << GetRowNum() << std::endl;
	std::cout << "ColNum: " << GetColNum() << std::endl;
	std::cout << "Rank: " << getRank() << std::endl;
	std::cout << "Trace: " << getTrace() << std::endl;
	if (GetRowNum() == GetColNum())
	{
		std::cout << "Square Matrix" << std::endl;
		std::cout << "Determinant: " << Determinant() << std::endl;
	}
	std::cout << "---------------------------------------" << std::endl;
	Print();
	std::cout << "***************************************" << std::endl;
}
template <typename T>
T Matrix<T>::Cofactor(uint32_t i, uint32_t j) const
{
	//求余子式
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::Cofactor: Only square matrix has cofactor!" << std::endl;
		return T();
	}
	if (i >= GetRowNum() || j >= GetColNum())
	{
		std::cout << "Error: Matrix::Cofactor: Index out of range!" << std::endl;
		return T();
	}
	Matrix<T> tmp(GetRowNum() - 1, GetColNum() - 1);
	for (uint32_t m = 0; m < GetRowNum(); m++)
	{
		for (uint32_t n = 0; n < GetColNum(); n++)
		{
			if (m < i && n < j)
			{
				tmp(m, n) = data[m][n];
			}
			else if (m < i && n > j)
			{
				tmp(m, n - 1) = data[m][n];
			}
			else if (m > i && n < j)
			{
				tmp(m - 1, n) = data[m][n];
			}
			else if (m > i && n > j)
			{
				tmp(m - 1, n - 1) = data[m][n];
			}
		}
	}
	return tmp.Determinant();
}
template <typename T>
Matrix<T>& Matrix<T>::Cofactorize()
{
	//求伴随矩阵
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::Adjoint: Only square matrix has adjoint matrix!" << std::endl;
		return *this;
	}
	Matrix<T> tmp(GetRowNum(), GetColNum());
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			tmp(i, j) = Cofactor(i, j) * ((i + j) % 2 == 0 ? (T)1 : (T)-1);

		}
	}
	*this = tmp;
	return *this;
}
template <typename T>
Matrix<T> Matrix<T>::GetCofactorization() const
{
	//求伴随矩阵
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::GetAdjoint: Only square matrix has adjoint matrix!" << std::endl;
		return Matrix<T>();
	}
	Matrix<T> tmp = *this;
	tmp.Cofactorize();
	return tmp;
}
template <typename T>
Matrix<T>& Matrix<T>::Adjoint()
{
	Cofactorize();
	Transpose();
	return *this;
}
template <typename T>
Matrix<T> Matrix<T>::GetAdjoint() const
{
	Matrix<T> tmp = *this;
	tmp.Cofactorize();
	tmp.Transpose();
	return tmp;
}