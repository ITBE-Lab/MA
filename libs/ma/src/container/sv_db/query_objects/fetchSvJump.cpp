#include "db_config.h"
#include "container/sv_db/query_objects/fetchSvJump.h"

using namespace libMA;

#ifdef WITH_PYTHON

using DBCon = SQLDB<MySQLConDB>;

void exportSvJump(py::module& rxPyModuleId)
{
	// export the SortedSvJumpFromSql class
	py::class_<SortedSvJumpFromSql<DBCon>>(rxPyModuleId, "SortedSvJumpFromSql")
		.def(py::init<const ParameterSetManager&, std::shared_ptr<SV_Schema<DBCon>>, int64_t>())
		.def(py::init<const ParameterSetManager&, std::shared_ptr<SV_Schema<DBCon>>, int64_t, int64_t, int64_t, uint32_t,
			uint32_t>())
		.def("has_next_start", &SortedSvJumpFromSql<DBCon>::hasNextStart)
		.def("has_next_end", &SortedSvJumpFromSql<DBCon>::hasNextEnd)
		.def("next_start_is_smaller", &SortedSvJumpFromSql<DBCon>::nextStartIsSmaller)
		.def("get_next_start", &SortedSvJumpFromSql<DBCon>::getNextStart)
		.def("get_next_end", &SortedSvJumpFromSql<DBCon>::getNextEnd);
} // function

#endif // WITH_PYTHON
