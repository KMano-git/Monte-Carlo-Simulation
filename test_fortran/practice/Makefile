# 1. 変数の定義
FC = gfortran
FFLAGS = -O3 -Wall
TARGET = simulation.out
SRCS = a.f90
OBJS = $(SRCS:.f90=.o)

# 2. 最終成果物の生成ルール（make と打つとここが実行される）
$(TARGET): $(OBJS)
	$(FC) $(OBJS) -o $(TARGET)

# 3. 各 .o ファイルの生成ルール
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# 4. 掃除用コマンド（make clean と打つと実行される）
clean:
	rm -f $(OBJS) $(TARGET) *.mod