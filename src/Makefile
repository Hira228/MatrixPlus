GCC=gcc
SRC=s21_matrix.cc
OBJ=s21_matrix.o
CFLAGS=--std=c++17 -lstdc++ -lm
TESTFLAGS=-lgtest -lgcov
GCOVFLAGS=--coverage
HTML=lcov -t test -o rep.info -c -d ./
OS = $(shell uname)

all: clean gcov_report

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt *.dSYM

test: s21_matrix_oop.a
	$(GCC) -g s21_matrix_oop_test.cc s21_matrix_oop.a $(CFLAGS) $(TESTFLAGS) -o test
	./test

s21_matrix_oop.a: clean
	$(GCC) $(GCOVFLAGS) -c $(SRC)
	ar rcs s21_matrix_oop.a $(OBJ)
	ranlib s21_matrix_oop.a


gcov_report: test
	$(HTML)
	genhtml -o report rep.info
	open report/index.html



style: 
	clang-format -i *.cpp *.h $(SRC_TEST)

$(CLEAN) :
	rm -rf tests *.a *.o $(TESTS).*o ./html_report/ gcov_report *.gcno *.gcda *.info *.out Build