SET CLIENT_ENCODING TO UTF8;
BEGIN;
CREATE TABLE network.LA_links
AS (
select
    (select distinct nodes.gid from network.LA_nodes nodes where nodes.geom = ST_StartPoint(ST_LineMerge(motorway.geom))) as fromid,
    (select distinct nodes.gid from network.LA_nodes nodes where nodes.geom = ST_EndPoint(ST_LineMerge(motorway.geom))) as toid,
    ST_length(ST_transform(motorway.geom, 26986)) as length
from network.LA_motorway motorway
union all
select
    (select distinct nodes.gid from network.LA_nodes nodes where nodes.geom = ST_EndPoint(ST_LineMerge(motorway.geom))) as fromid,
    (select distinct nodes.gid from network.LA_nodes nodes where nodes.geom = ST_StartPoint(ST_LineMerge(motorway.geom))) as toid,
    ST_length(ST_transform(motorway.geom, 26986)) as length
from network.LA_motorway motorway
);
ALTER TABLE network.LA_links ADD COLUMN type character varying(15) DEFAULT 'HIGHWAY';
ALTER TABLE network.LA_links ADD COLUMN capacity double precision DEFAULT 8000.0;
ALTER TABLE network.LA_links ADD COLUMN freespeed double precision DEFAULT 30.0;
ALTER TABLE network.LA_links ADD COLUMN gid SERIAL NOT NULL;
ALTER TABLE network.LA_links ADD COLUMN id SERIAL NOT NULL;
COMMIT;
