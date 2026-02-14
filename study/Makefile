# C++コンパイル用Makefile
# 使い方: make ファイル名（拡張子なし）
# 例: make ex1  →  g++ -std=c++17 -O2 -o ex1 ex1.cpp

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

# パターンルール: % は任意のファイル名にマッチ
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# クリーンアップ（実行ファイルを削除）
clean:
	find . -type f -executable -name '[!.]*' -delete

.PHONY: clean
