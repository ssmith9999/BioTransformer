#/usr/bin/bash

root_dir=$(readlink -f "$(dirname $0)")
echo $root_dir 

#Cleaning
echo "Cleaning..."
mvn clean

#Adding phaseIIFilter
echo -e "\n\nAdding Phase II Filter JAR..."
mvn org.apache.maven.plugins:maven-install-plugin:3.0.0-M1:install-file -Dfile="${root_dir}/lib/phaseIIFilter-1.0.1.jar" -DgroupId=djoumbou -DartifactId=phaseIIFilter -Dversion=1.0.1 -Dpackaging=jar


#Adding Cypreact
echo -e "\n\nAdding CypReact JAR..."
mvn org.apache.maven.plugins:maven-install-plugin:3.0.0-M1:install-file -Dfile="${root_dir}/lib/cypreact.jar" -DgroupId=wishartlab -DartifactId=cypreact -Dversion=1.0.0-SNAPSHOT -Dpackaging=jar

#Building package
echo -e "\n\nBuilding Package..." 
mvn package

#Moving jar file
echo -e "Moving .jar file to bin/ ..."
mv target/*.jar bin/

echo "Done!"


