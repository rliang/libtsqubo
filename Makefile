CFLAGS ?= -Wall -Wextra -Werror -Og -g

default: test_nopq test_pq

test_%: test-%.out
	curl -fsSL https://web.stanford.edu/~yyye/yyye/Gset/G1 | ./$^

test-nopq.out: test.c tsqubo.h
	$(CC) $(CFLAGS) test.c -o $@

test-pq.out: test.c tsqubo.h
	$(CC) $(CFLAGS) -DTSQUBO_SPARSE test.c -o $@

.INTERMEDIATE: test-nopq.out test-pq.out
