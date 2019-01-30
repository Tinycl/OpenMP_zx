#include <windows.h>
#include <iostream>
#include <omp.h>

using namespace std;

void SimpleWithOMP(int n, float *a, float *b)
{
	int i;
	int iam, nthd, iend, istart;
	#pragma omp parallel for //default(shared) private(iam,nthd, istart, iend)
	for (i = 1; i < n; i++)
	{
		b[i] = (a[i] + a[i - 1]) / 2.0 * 3.0;
		//printf("threadid is %d  i is %d \n",omp_get_thread_num(),i); //��ӡcore number
	}
}
void SimpleNoOMP( int n, float *a, float *b)
{
	int i;
	for (i = 1; i < n; i++)
	{
		b[i] = (a[i] + a[i - 1]) / 2.0 * 3.0;
	}
}

void subdomain(float *x, int istart, int ipoints)
{
	int i;
	for (i = 0; i < ipoints; i++)
	{
		x[istart + i] = 123.456;
		printf("threadid is %d \n",omp_get_thread_num());
	}
}
//#pragma omp parallel ʵ�� for ��Ĳ���
void sub(float *x, int npoints)
{
	int iam, nt, ipoints, istart;
	#pragma omp parallel default(shared) private(iam,nt, ipoints, istart)
	{
		iam = omp_get_thread_num(); //�ĸ�coreִ��  ������Ϊcore���Ϳ��Բ����Ż����߳���һ��
		nt = omp_get_num_threads(); //�ܵ�core����
		ipoints = npoints / nt; //size of partition �������
		istart = iam * ipoints; //starting arrray index
		if (iam == nt - 1)   //last thread may do more or less
		{
			ipoints = npoints - istart;
		}
		subdomain(x,istart,ipoints);
	}
}
//parallel ʵ�� ���� for �������� �ȼ��� parallel for
void array_iterate_cal(float *array, int length)
{
	int iam, istart, iend, nthrd, iter;
	#pragma omp parallel default(shared) private(iam,istart, iend, nthrd, iter) num_threads(4)
	{
		iam = omp_get_thread_num();  //�������� ����ִ�е��߳�
		nthrd = omp_get_num_threads(); //�������� �����ٸ��߳�
		iend = length / nthrd;   //������Ԫ�� �ָ�ÿ���̴߳���ÿһ�εĳ���
		istart = iam * iend;    //ÿһ�εĿ�ʼ
		if (iam == nthrd - 1)  //���һ���̣߳������Ԫ�ؿ����� С��iend
		{
			iend = length - istart;
		}
		for (iter = 0; iter < iend; iter++)
		{
			array[istart + iter] = 1.0*2.0 + 1.0;
			printf("threadid is %d  iter is %d \n", omp_get_thread_num(), istart + iter);
		}
	}
}

void work1()
{

}
void work2()
{

}
static omp_lock_t g_lock;
int main()
{
	int arraysize = 10;// 2592 * 1944;// 2592 * 1944; //2048*1536
	float *a = new float[arraysize];
	for (int i = 0; i < arraysize; i++)
	{
		a[0] = i * 1.0;
	}
	float *b = new float[arraysize];
#ifndef _OPENMP
	fprintf(stderr, "OpenMP not supported");
#endif
	omp_set_num_threads(4); //����4��core���д��� 4core cpu
	//�����Ż�����
	LARGE_INTEGER timestart;
	LARGE_INTEGER timeend;
	LARGE_INTEGER frequency;
	double quadpart = 0;
	double elapsed = 0;
	printf("\/************optimation difference************\/\n");
	
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&timestart);
	quadpart = (double)frequency.QuadPart; 
	SimpleWithOMP(arraysize,a,b);
	QueryPerformanceCounter(&timeend);
	elapsed = (timeend.QuadPart - timestart.QuadPart) / quadpart; 
	printf("omp optimation for loop NO : %f s \n", elapsed);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&timestart);
	quadpart = (double)frequency.QuadPart;
	SimpleNoOMP(arraysize, a, b);
	QueryPerformanceCounter(&timeend);
	elapsed = (timeend.QuadPart - timestart.QuadPart) / quadpart;
	printf("omp optimation for loop YES : %f s \n", elapsed);
	if (a != NULL) {
		delete[] a;
	}
	if (b != NULL) {
		delete[] b;
	}
	
	//������䣬������
	printf("\/************array partition************\/\n");
	float array[10];
	sub(array, 10);
	//Ƕ��
	printf("\/************nest test************\/\n");
	float *array1 = new float[2555*2555];
	double ompstarttime = 0, ompendtime = 0;
	omp_set_nested(1); //����һ�㲢��Ƕ��
	ompstarttime = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < 2555; i++) {
		//#pragma omp parallel for Ƕ��һ��ѭ��ʱ�� ֻ����������ָ����������Ż������У�
		for(int j = 0; j < 2555; j++){
			array1[i*2555+j] = 2.0*1.5 + 1;
		}
	}
	ompendtime = omp_get_wtime();
	printf("use 1 nest to cal 2 dim : %f s \n", ompendtime - ompstarttime);
	
	omp_set_nested(0); //�����ò���Ƕ��
	ompstarttime = omp_get_wtime();
	for (int i = 0; i < 2555; i++) {
		for (int j = 0; j < 2555; j++) {
			array1[i*2555+j] = 2.0*1.5 + 1;
		}
	}
	ompendtime = omp_get_wtime();
	printf("don't use nest to cal 2 dim : %f s \n", ompendtime - ompstarttime);
	delete[] array1;
	//���ö�̬����
	printf("\/************set dynamic************\/\n");
	float array2[10];
	omp_set_dynamic(0); //���ö�̬����
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < 10; i++) {
		array2[i] = 1.0*2.0 + 1.0;
		printf("threadid is %d \n", omp_get_thread_num());
	}

	//�߳�ͬ��֮omp_lock
	printf("\/************mutex lock************\/\n");
	omp_init_lock(&g_lock);  //��ʼ��lock
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < 10; i++) {
		omp_set_lock(&g_lock);  //��û�����
		printf("threadid is %d +\n", omp_get_thread_num()); //����3��printf������
		printf("threadid is %d -\n", omp_get_thread_num());
		printf("threadid is %d *\n", omp_get_thread_num());
		omp_unset_lock(&g_lock);  //�ͷŻ�����
	}
	omp_destroy_lock(&g_lock);

	//nowait
	printf("\/************nowait************\/\n");
	#pragma omp parallel 
	{
		//pragma omp for Ĭ�Ϻ����и���ʽBarrier
		#pragma omp for nowait     // ����for block���Բ���ִ�м� +-+-- �����ǵȵ���һ��for��ɺ�ִ����һ��for +++---
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d +\n", omp_get_thread_num()); 
		}
		#pragma omp for 
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num());
		}
	}
	
	//barrier
	printf("\/************barrier************\/\n");
	#pragma omp parallel 
	{
		
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d +\n", omp_get_thread_num());
		}
		#pragma omp barrier  //����ʾ����barrier�� ����for�鲢��ִ�С���ʾ����barrier����ִ�����һ��for����ִ�еڶ���for
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num());
		}
	}

	//master
	printf("\/************master************\/\n");
	#pragma omp parallel 
	{
		#pragma omp master  //��for���п�ֻ�����߳�ִ��
		{
			for (int i = 0; i < 10; i++) {
				printf("threadid is %d +\n", omp_get_thread_num()); //��ӡ10��
			}
		}
		
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num()); //��ӡ10*core��
		}
	}

	//pareallel -> parallel for
	printf("\/************parallel -> parallel for************\/\n");
	array_iterate_cal(array2,10);

	//sections
	printf("\/************sections************\/\n");
	#pragma omp parallel sections
	{
		//����section for �ǲ���ִ�е�
		#pragma omp section  //ֻ��һ���߳�printf10�� 
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d +\n", omp_get_thread_num());
		}
		#pragma omp section //ֻ������һ���߳�printf10��
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num());
		}
	}
	//atomic
	printf("\/************atomic************\/\n");
	int sum = 0;
	#pragma omp parallel for
	for (int i = 1; i <= 50; i++)
	{
		//#pragma omp atomic
		sum += i;
		printf("threadid is %d \n", omp_get_thread_num());
	}
	
	//single
	printf("\/************single************\/\n");
	#pragma omp parallel
	{
		#pragma omp single  //���д������ڣ�ֻ��һ���̣߳�����ִ�иô���
		printf("Beginning work1 \n");
		work1();
		#pragma omp single
		printf("Finishing work1 \n");

		#pragma omp single nowait
		printf("Finish work1 and beginning work2 \n");
		work2();
	}
	system("pause");
	return 0;
}