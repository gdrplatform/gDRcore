# sudo apt-get install plantuml
# pretty old version
#jar_path="/usr/share/platinuml/platinuml.jar"
# https://github.com/plantuml/plantuml/releases/download/v1.2023.13/plantuml-1.2023.13.jar
jar_path="./plantuml-1.2023.13.jar"
#java -Djava.awt.headless=true -jar $jar_path *.puml -tsvg -o ../images/
java -Djava.awt.headless=true -jar $jar_path *.puml -tpng -o ../images/
