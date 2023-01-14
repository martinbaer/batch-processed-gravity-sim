if compgen -G "objects/*.o" > /dev/null; then \
	mv objects/*.o .; \
fi
if compgen -G "objects/*.o" > /dev/null; then \
	rm objects/*.o; \
fi