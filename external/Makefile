all: catch nlohmann_json cereal Clara

.PHONY: all

nlohmann_json: nlohmann_json-3.1.2.tar
	tar -mxf $<
	ln -s nlohmann_json-*/ $@

Clara: Clara-ba5485c.tar
	tar -mxf $<
	ln -s Clara-*/ $@
