CREATE TABLE Sequence(
	ID INTEGER PRIMARY KEY AUTOINCREMENT, --Auto incrementing key value
	Plate INT,
	Clone TEXT,
	Primer TEXT check(Primer='F' or Primer='R'),
	Run_date TEXT,
	Student TEXT,
	Instructor TEXT,
	Institution TEXT,
	Quality INT check(Quality=0 or Quality=1),
	FASTA TEXT NOT NULL,
	ABI TEXT NOT NULL,
	Pursue INT check(Pursue=0 or Pursue=1 or Pursue is NULL),
	'Comment' TEXT
);

CREATE TABLE Blast(
	ID INT check(ID is NOT NULL),
	Genome TEXT,
	Organism TEXT,
	Query_from INT,
	Query_to INT,
	Subject_from INT,
	Subject_to INT,
	Identity INT,
	Pursue INT,
	Similarity TEXT,
	FOREIGN KEY(ID) REFERENCES Sequence(ID)
);

CREATE TABLE toBLAST(
	ID INT,
	FASTA TEXT  NOT NULL,
	FOREIGN KEY(ID) REFERENCES Sequence(ID)
);

-- Every new high quality insert to the Sequence table will also have a record in the toBLAST table
CREATE TRIGGER addBlast AFTER INSERT ON Sequence WHEN new.Quality=1
BEGIN
	INSERT INTO toBlast VALUES(new.ID, new.FASTA);
END;

-- Every Blast parsing result will be updated on the Sequence table whether to pursue or not
CREATE TRIGGER not_pursue AFTER INSERT ON Blast WHEN new.Pursue=0
BEGIN
	UPDATE Sequence SET Pursue=0 WHERE Sequence.ID=new.ID;
END;

CREATE TRIGGER pursue AFTER INSERT ON Blast WHEN new.Pursue=1
BEGIN
	UPDATE Sequence SET Pursue=1 WHERE Sequence.ID=new.ID;
ENd;

CREATE VIEW Rerun AS
	select DISTINCT L.Plate, L.Clone, L.Primer
	from
		(select DISTINCT Plate, Clone, Primer
		from Sequence
		where Sequence.Quality=1) as S,
		(select  DISTINCT Plate, Clone, Primer, Run_date
		from Sequence
		where Quality=0) as L
	where not (S.Plate = L.Plate and S.Clone=L.Clone and S.Primer=L.Primer);

CREATE VIEW Pursue AS
	select *
	from Sequence
	where Pursue = 1;
