default: test.out G1
	./test.out < G1

test.out: tsqubo.h test.c Makefile
	$(CC) -Wall -Wextra -Werror -Og -g test.c -o $@

G1:
	curl -fLO https://web.stanford.edu/~yyye/yyye/Gset/G1