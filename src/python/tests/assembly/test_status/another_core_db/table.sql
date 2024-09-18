CREATE TABLE meta (
  meta_id    INTEGER NOT NULL,
  species_id INTEGER DEFAULT 1,
  meta_key   VARCHAR(64) NOT NULL,
  meta_value VARCHAR(255) DEFAULT '',
  PRIMARY KEY (species_id, meta_key, meta_value)
);
CREATE INDEX meta_id_idx ON meta (meta_id);