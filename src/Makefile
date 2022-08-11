CC = g++
CPPFLAGS = -lgtest -std=c++17
CHECKFLAGS = -Wall -Werror -Wextra -g

all: clean s21_matrix_oop.a test

test:
	@${CC} test.cc s21_matrix_oop.cc ${CHECKFLAGS} ${CPPFLAGS}
	@./a.out

s21_matrix_oop.a:
	@$(CC) $(CHECKFLAGS) s21_matrix_oop.cc -c
	@ar r s21_matrix_oop.a s21_matrix_oop.o

clean:
	@/bin/rm -rf *.o *.a matrix *.gcno *gcda report *.info matrix *.out *.dSYM *.cfg

gcov_report: test.cc s21_matrix_oop.cc
	@$(CC) --coverage -lcheck $^ ${CPPFLAGS} ${CHECKFLAGS}
	@./a.out
	@lcov -t a.out -o rep.info --no-external -c -d .
	@genhtml -o report rep.info
	@open report/index.html

check:
	@cppcheck --enable=all --suppress=missingIncludeSystem . *.cc
	@cp ../materials/linters/CPPLINT.cfg .
	@python3 ../materials/linters/cpplint.py --extensions=cc *.cc *.h
	@rm -rf CPPLINT.cfg
	@make test
	@CK_FORK=no leaks --atExit -- ./a.out

git:
	git add .
	git commit -m 'a'
	git push origin develop
