DB_PATH = "Sequence Files/Seq.db"

# Quality threshold for checking sequecing quality
MIN_NUCLEOTIDES = 300
GOOD_THRESHOLD = 30
UNCERTAIN_THRESHOLD = 20

# Filter out trace files not non-samples
BLACK_LIST = ['water', 'WATER', 'STANDARD', 'BigDye', 'BUFFER', 'DNA1F', 'DNA1R', 'DNA2F', 'DNA2R', 'DNA3F', 'DNA3R', 'DNAA7F', 'DNAA7R', 'DNAA8F', 'DNAA8R', 'DNAA9F', 'DNAA9R', 'control1', 'control2', 'control3', 'control4', 'NBCCONTROL', 'CONTROL', 'empty', 'unknown', 'Test']

# Parameters for the blast search
E_THRESHOLD = 1e-50
alignments = 100

# List of super-colonizer and non-colonizer for the blast to distinguish pursue or not
COLONIZER = ['brassicacearum', 'NFM421', 'CP002585', 'F113', 'CP003150', 'Q8r1-96']
NON_COLONIZER = ['Pf-5', 'CP000076', 'Pf0-1', 'CP000094', 'SBW25', 'AM181176', 'Q287']

# Database schema
schema = """
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
    E_value REAL,
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

--Delete from toBlast table when an entry is removed from the Sequence table
CREATE TRIGGER remove_pursue AFTER DELETE ON Sequence
BEGIN
    DELETE FROM toBlast WHERE toBlast.ID=old.ID;
END;

-- Every Blast parsing result will be updated on the Sequence table whether to pursue or not
CREATE TRIGGER update_not_pursue AFTER INSERT ON Blast WHEN new.Pursue=0
BEGIN
    UPDATE Sequence SET Pursue=0 WHERE Sequence.ID=new.ID;
END;

CREATE TRIGGER update_pursue AFTER INSERT ON Blast WHEN new.Pursue=1
BEGIN
    UPDATE Sequence SET Pursue=1 WHERE Sequence.ID=new.ID;
END;

--Delete the record from toBlast table if the Blast result has been inserted
CREATE TRIGGER toBlast_clean_up AFTER INSERT ON Blast
BEGIN
    DELETE FROM toBlast WHERE toBlast.ID=new.ID;
END;

CREATE VIEW LowQuality AS
    select DISTINCT Plate, Clone, Primer
    from Sequence
    where Quality = 0;

CREATE VIEW HighQuality As
    select DISTINCT Plate, Clone, Primer
    from Sequence
    where Quality = 1;

--Cheap trick where the entries are concatenated and compared to each others
CREATE VIEW Rerun AS
    select Plate, Clone, Primer
    from LowQuality
    where Plate||Clone||Primer not in (select Plate||Clone||Primer from HighQuality);

CREATE VIEW Pursue AS
    select *
    from Sequence
    where Pursue = 1;
"""

# Database insert commands

SQL_COMMAND = {
    'Sequence': 'INSERT INTO Sequence(Plate, Clone, Primer, Run_date, Student, Instructor, Institution, Quality, FASTA, ABI, "Comment") VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)',
    'Blast': 'INSERT INTO Blast(ID, Genome, Organism, E_value, Query_from, Query_to, Subject_from, Subject_to, Identity, Pursue, Similarity) VALUES(?,?,?,?,?,?,?,?,?,?,null)',
    'remove_last': 'DELETE FROM Sequence WHERE ID = (SELECT MAX(ID) FROM Sequence)'}
