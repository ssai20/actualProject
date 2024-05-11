DROP DATABASE IF EXISTS solution_base;
CREATE DATABASE solution_base;
USE solution_base;

CREATE TABLE solution (
id INT(11) NOT NULL AUTO_INCREMENT,
title VARCHAR(50) NOT NULL,
PRIMARY KEY (id)
);

INSERT INTO solution value(NULL, "Таблица 1");


SELECT title FROM solution;