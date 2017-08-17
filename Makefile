.PHONY: clean All

All:
	@echo "----------Building project:[ ANITA - Debug ]----------"
	@"$(MAKE)" -f  "ANITA.mk"
clean:
	@echo "----------Cleaning project:[ ANITA - Debug ]----------"
	@"$(MAKE)" -f  "ANITA.mk" clean
