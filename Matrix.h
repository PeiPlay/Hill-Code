#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <initializer_list>



//�����࣬���ڴ洢�������ݺͽ��о�������
//����Ϊģ���࣬���Դ洢�������͵�����

template <typename T>
class Matrix
{
public:
	Matrix(uint32_t m, uint32_t n);									//���캯����mΪ������nΪ����
	Matrix(const Matrix<T>& mat);									//�������캯��
	Matrix(const T* mat, uint32_t m, uint32_t n);					//���캯����matΪ����ָ�룬mΪ������nΪ����
	Matrix(std::initializer_list<std::initializer_list<T>> list);	//�����ų�ʼ������

	//���������
	Matrix<T>& operator=(const Matrix<T>& mat);			//��ֵ���������
	Matrix<T> operator+(const Matrix<T>& mat) const;	//�ӷ����������
	Matrix<T> operator-(const Matrix<T>& mat) const;	//�������������
	Matrix<T> operator*(const Matrix<T>& mat) const;	//�˷����������
	Matrix<T> operator*(const T& num) const;			//�������������
	Matrix<T>& operator+=(const Matrix<T>& mat);		//�ӷ���ֵ���������
	Matrix<T>& operator-=(const Matrix<T>& mat);		//������ֵ���������
	Matrix<T>& operator*=(const Matrix<T>& mat);		//�˷���ֵ���������
	Matrix<T>& operator*=(const T& num);				//���˸�ֵ���������
	bool operator==(const Matrix<T>& mat) const;		//������������
	bool operator!=(const Matrix<T>& mat) const;		//�������������
	

	//�����������
	void Add(const Matrix<T>& mat);		//����ӷ�
	void Sub(const Matrix<T>& mat);		//�������
	void Mul(const Matrix<T>& mat);		//����˷�
	void Mul(const T& num);				//��������

	void setUnit(uint32_t m, uint32_t n, T u)
	{
		data[m][n] = u;
	}

	T getUnit(uint32_t m, uint32_t n) const
	{
		return data[m][n];
	}

	friend void GetRandomInvertibleMatrix(uint32_t seed)	//��������������
	{

	}

	void RowExchange(uint32_t i, uint32_t j);				//�����н���
	void ColExchange(uint32_t i, uint32_t j);				//�����н���
	void RowAdd(uint32_t i, uint32_t j, T k);				//�����б���
	void ColAdd(uint32_t i, uint32_t j, T k);				//�����б���
	

	//�������ݷ���
	T& operator()(uint32_t i, uint32_t j);					//���ʾ����i�е�j�е�����
	const T& operator()(uint32_t i, uint32_t j) const;		//���ʾ����i�е�j�е�����
	Matrix<T> GetRow(uint32_t i) const;						//��ȡ�����i��
	Matrix<T> GetCol(uint32_t j) const;						//��ȡ�����j��
	uint32_t getRank() const;								//��ȡ������
	uint32_t getNullity() const;							//��ȡ�������
	T getTrace() const;										//��ȡ����
	T getNorm() const;										//��ȡ������


	//������ز���
	Matrix<T> GetEchelonForm() const;			//��ȡ���������ʽ
	Matrix<T>& ToEchelonForm();					//����Ϊ�����ʽ
	Matrix<T> GetReducedEchelonForm() const;	//��ȡ�����������ʽ
	Matrix<T>& ToReducedEchelonForm();			//����Ϊ�������ʽ


	Matrix<T>& Transpose();						//����ת��
	Matrix<T> GetTranspose() const;				//��ȡ����ת��

	Matrix<T>& Inverse();						//��������
	Matrix<T> GetInverse() const;				//��ȡ���������
	bool isInvertible() const;					//�жϾ����Ƿ����

	T Cofactor(uint32_t i, uint32_t j) const;	//����������ʽ

	Matrix<T>& Cofactorize();					//����������ʽ����
	Matrix<T> GetCofactorization() const;		//��ȡ��������ʽ����

	Matrix<T>& Adjoint();						//������������
	Matrix<T> GetAdjoint() const;				//��ȡ����������

	T Determinant() const;						//����������ʽ

	//���������Ϣ
	void Print() const;							//��ӡ����
	void ShowInfo() const;						//��ʾ������Ϣ

	uint32_t GetRowNum() const;					//��ȡ��������
	uint32_t GetColNum() const;					//��ȡ��������


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

	//��������
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
	//�ڲ�����д�����

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
void Matrix<T>::RowExchange(uint32_t i, uint32_t j)				//�����н���
{
	//������i����j
	std::vector<T> temp = data[i];
	data[i] = data[j];
	data[j] = temp;
}
template <typename T>
void Matrix<T>::ColExchange(uint32_t i, uint32_t j)				//�����н���
{
	//������i����j
	for (uint32_t k = 0; k < data.size(); k++)
	{
		T temp = data[k][i];
		data[k][i] = data[k][j];
		data[k][j] = temp;
	}
}
template <typename T>
void Matrix<T>::RowAdd(uint32_t i, uint32_t j, T k)				//�����б���
{
	//����j��k���ӵ���i��
	for (uint32_t l = 0; l < data[i].size(); l++)
	{
		data[i][l] += data[j][l] * k;
	}
}
template <typename T>
void Matrix<T>::ColAdd(uint32_t i, uint32_t j, T k)				//�����б���
{
	//����j��k���ӵ���i��
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
		//�ҵ���j�е�i�м���֮������Ԫ�����ڵ���
		uint32_t maxRow = i;
		for (uint32_t k = i + 1; k < m; k++)
		{

			if (Abs(data[k][j]) > Abs(data[maxRow][j]))
			{
				maxRow = k;
			}
		}
		//������Ԫ��Ϊ0������jȫΪ0��������һ��
		if (data[maxRow][j] == 0)
		{
			j++;
		}
		else
		{
			//������i�к͵�maxRow��
			if (i != maxRow)
			{
				for (uint32_t k = j; k < n; k++)
				{
					std::swap(data[i][k], data[maxRow][k]);
				}
			}
			//����i�е���Ԫ������Ϊ1
			T scale = data[i][j];
			for (uint32_t k = j; k < n; k++)
			{
				data[i][k] /= scale;
			}
			//����i�����µĵ�j��Ԫ����Ϊ0
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
		//�ҵ���j�е�i�м���֮������Ԫ�����ڵ���
		uint32_t maxRow = i;
		for (uint32_t k = i + 1; k < m; k++)
		{
			if (Abs(data[k][j]) > Abs(data[maxRow][j]))
			{
				maxRow = k;
			}
		}
		//������Ԫ��Ϊ0������jȫΪ0��������һ��
		if (data[maxRow][j] == 0)
		{
			j++;
		}
		else
		{
			//������i�к͵�maxRow��
			if (i != maxRow)
			{
				for (uint32_t k = j; k < n; k++)
				{
					std::swap(data[i][k], data[maxRow][k]);
				}
			}
			//����i�е���Ԫ������Ϊ1
			T scale = data[i][j];
			for (uint32_t k = j; k < n; k++)
			{
				data[i][k] /= scale;
			}
			//����i�����µĵ�j��Ԫ����Ϊ0
			for (uint32_t k = i + 1; k < m; k++)
			{
				T scale = data[k][j];
				for (uint32_t l = j; l < n; l++)
				{
					data[k][l] -= scale * data[i][l];
				}
			}
			//����i�����ϵĵ�j��Ԫ����Ϊ0
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
	//ʵ�־���ת��
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
	//ʵ�־���������ʽ
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
	//ʵ���жϾ����Ƿ����
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
	//ʵ�־�������
	if (GetRowNum() != GetColNum())
	{
		std::cout << "Error: Matrix::Inverse: Only square matrix can be inverted!" << std::endl;
		return *this;
	}
	//ʹ�ø�˹��Ԫ������
	Matrix<T> augmentation(GetRowNum(), GetColNum() * 2);	//�������
	//��ʼ���������
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		for (uint32_t j = 0; j < GetColNum(); j++)
		{
			augmentation(i, j) = data[i][j];
		}
		augmentation(i, GetColNum() + i) = 1;
	}
	//��˹��Ԫ
	for (uint32_t i = 0; i < GetRowNum(); i++)
	{
		//�ҵ���i�е���Ԫ
		uint32_t max_row = i;
		for (uint32_t j = i + 1; j < GetRowNum(); j++)
		{
			if (augmentation(j, i) > augmentation(max_row, i))
			{
				max_row = j;
			}
		}
		//������i�к͵�max_row��
		for (uint32_t j = 0; j < GetColNum() * 2; j++)
		{
			std::swap(augmentation(i, j), augmentation(max_row, j));
		}
		//����i�е���Ԫ����Ϊ1
		T scale = augmentation(i, i);
		for (uint32_t j = 0; j < GetColNum() * 2; j++)
		{
			augmentation(i, j) /= scale;
		}
		//����i�е�����Ԫ����Ϊ0
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
	//�����������Ұ벿����Ϊ���
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
	//ʵ�־�������
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
	//��ӡ����
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
	//��ӡ������Ϣ
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
	//������ʽ
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
	//��������
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
	//��������
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