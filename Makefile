
doc:
	R -e 'devtools::document()'

test:
	R -e 'devtools::test()'

build:
	R -e 'devtools::build()'

install:
	R -e 'devtools::document()'
	R -e 'devtools::install()'

check:
	R -e 'devtools::test()'
