/* We keep the configuration in a property tree.
 */
#include "configure.h"
#include "global.h" // due to xGlobalConfiguration

/* Platform and OS informing macros.
 * For the 2 directory related functions.
 */
#include <boost/predef.h>

// #define assemblyFolder "../../../Data/database/Assembly/All/"//should be: "../../../Data/database/Assembly"
// #define assemblyFolder_NCBI "genomes/ASSEMBLY_REPORTS/All/"
// #define taxDump_Folder_local "../../../Data/database/taxdmp"

/* The standard configuration represented as XML and parsed during startup.
 */
std::string sGlobalStandardConfiguration (
"<configure>"
	"<log>"	// configuration for logging
		"<trivial>" // boost trivial logging
			"<trace>false</trace>" 
			"<debug>false</debug>"
			"<info>true</info>"
			"<warning>true</warning>"
			"<error>true</error>"
			"<fatal>true</fatal>"
		"</trivial>"
	"</log>"
R"(	<global>
		<prefix>)" 
#if _MSC_VER
			"F:\\BioSolution"
#else
			"/opt/bio" // "/usr/home/itbe"
#endif
		R"(</prefix>
		<temporary>
			<folder>tmp</folder>
		</temporary>
	</global>
	<NCBI>
		<prefix>NCBI</prefix>
		<server>
			<ftp>
				<hostname>ftp.ncbi.nlm.nih.gov</hostname>
			</ftp>
		</server>	
		<cache>
			<folder>cache</folder>
			<xml>
				<subfolder>xml</subfolder>
			</xml>
			<fasta>
				<subfolder>fasta</subfolder>
			</fasta>
		</cache>
		<assemblies>
			<folder>assemblies</folder>
			<sources>
				<location>genomes/ASSEMBLY_REPORTS/All</location> 
				<folder>sources</folder>
			</sources>
			<preassembled>
				<subfolder>preassembled</subfolder>
			</preassembled>
			<cache>)"
#if _MSC_VER
				"<folder>\\BioSolution\\NCBI\\assemblies</folder>"
#else
				"<folder>/home/markus/BioSolution/assemblies</folder>"//"<folder>/opt/bio/NCBI/assemblies</folder>" // "<folder>/usr/home/itbe/NCBI/assemblies</folder>"
#endif 
			R"(	<fm-indices>
					<subfolder>fm-indices</subfolder>
				</fm-indices>
				<packs>
					<subfolder>packs</subfolder>
				</packs>
			</cache>
		</assemblies>
		<database>
			<prefix>database</prefix>
			<name>NCBI.db</name>
			<taxonomy>
				<prefix>taxonomy</prefix>
			</taxonomy>
			<sources>
				<input> <location>pub/taxonomy/taxdump.tar.gz</location> <name>taxdump.tar.gz</name> </input>
				<input> <location>/gene/DATA/gene2accession.gz</location> <name>gene2accession</name> </input>
				<input> <location>/gene/DATA/gene_info.gz</location> <name>gene_info</name> </input>
				<input> <location>/gene/DATA/gene_neighbors.gz</location> <name>gene_neighbors</name> </input>
				<input> <location>/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt</location> <name>assembly_summary_refseq.txt</name> </input>
				<input> <location>/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt</location> <name>assembly_summary_genbank.txt</name> </input>
			</sources>
		</database>
	</NCBI>
)"
"</configure>"
); // configuration string

/* The trivial log filter function.
 * Currently implemented as global function without class scope.
 */
bool trivialLogFilter( boost::log::value_ref< boost::log::trivial::severity_level, boost::log::trivial::tag::severity > const& rxSeveritylevel, // severity of message
					   decltype( Configuration::aLogLevelSwitches ) const &aLogLevelSwitches
				     )
{
	auto eSeveritylevel = rxSeveritylevel.get(); // extract the enumeration value from the value reference 	
	assert( eSeveritylevel <= boost::log::trivial::fatal ); // during debugging check that the index is fine.
	return aLogLevelSwitches[eSeveritylevel]; // deliver the status for the requested logging level.
} // lambda

/* Initializes the behavior of trivial logging.
 */
void Configuration::vInitializeLogging( void )
{
	/* Initialize the array that comprises the filter switches.
	 */
	aLogLevelSwitches[boost::log::trivial::trace] = xPropertyTree.get<bool>( "configure.log.trivial.trace" );
	aLogLevelSwitches[boost::log::trivial::debug] = xPropertyTree.get<bool>( "configure.log.trivial.debug" );
	aLogLevelSwitches[boost::log::trivial::info] = xPropertyTree.get<bool>( "configure.log.trivial.info" );
	aLogLevelSwitches[boost::log::trivial::warning] = xPropertyTree.get<bool>( "configure.log.trivial.warning" );
	aLogLevelSwitches[boost::log::trivial::error] = xPropertyTree.get<bool>( "configure.log.trivial.error" );
	aLogLevelSwitches[boost::log::trivial::fatal] = xPropertyTree.get<bool>( "configure.log.trivial.fatal" );
	
	/* Set the handler filter function that is called in the context of logging. 
	 */
	boost::log::core::get()->set_filter
	(	/* Taken from the boost log examples. Works only with boost::phoenix::bind and not with std::bind, for whaterver reasons.
		 * Extended form of a simple filer like: boost::log::trivial::severity >= boost::log::trivial::info
		 */
		boost::phoenix::bind
		(
			&trivialLogFilter, // filter as lambda does not work (presumably because lambdas are "function-typeless" functors in C++)
			boost::log::trivial::severity.or_none(), // severity of message
			this->aLogLevelSwitches	// array with switch that decide the logging behavior
		) // function call (phoenix::bind)
	); // function call (set_filter)
} // method

/* The global configuration object, which is initialized at program start.
 * Maybe there is better location for for this function ...
 */
Configuration xGlobalConfiguration( sGlobalStandardConfiguration );

/* Delivers the home directory of the current user for Linux and Windows.
 */
boost::filesystem::path sGetHomeDirectory() 
{
#if BOOST_OS_LINUX == 1
	std::string sHomePath = getenv( "HOME" ) != NULL ? getenv( "HOME" ) : "";
#elif BOOST_OS_WINDOWS == 1
	std::string sHomePath = getenv( "USERPROFILE" ) != NULL 
		? getenv( "USERPROFILE" )
		: (getenv( "HOMEDRIVE" ) != NULL && getenv( "HOMEPATH" ) != NULL)
			? std::string( getenv( "HOMEDRIVE" ) ).append( getenv( "HOMEPATH" ) )
			: "";
#else
	#error Get home directory undefined for current OS. 
#endif
	
	return boost::filesystem::path( sHomePath );
} // function

/* Delivers (creates, if necessary) a empty directory on the foundation of rsPath.
 * If the requested directory exists already, we create a fresh one by following a numbering scheme.
 * Maybe there is better location for for this function ...
 */
boost::filesystem::path sMakeDirectoryIfNotExists( const boost::filesystem::path &rsPath )
{
	boost::filesystem::path sSearchPath( rsPath );
	size_t uiCounter = 1;
	
	while( true )
	{	/* Check whether an directory with the requested path already exists
		 */
		if ( boost::filesystem::exists( sSearchPath ) )
		{
			if ( boost::filesystem::is_empty( sSearchPath ) )
			{	/* Directory exists and is empty.
				 */
				break;
			} // if
			else
			{	/* Directory exists and is NOT empty. Create new search path. (Waring += is right associative!)
				 */
				sSearchPath = (( boost::filesystem::path( rsPath ) += "_" ) += std::to_string( uiCounter++ ));
			} // else
		} // if
		else
		{	/* There is no directory with path rsPath.
			 * Create such a directory.
			 */
			boost::system::error_code xErrorCode;
			if (! boost::filesystem::create_directories(  sSearchPath,  xErrorCode ) )
			{
				throw std::runtime_error( std::string( "Can't create directory " + sSearchPath.string() + " due to code " + xErrorCode.message() ).c_str() );
			} // if
			
			/* Directory didn't exist, but successfully created it.
			 */
			break;
		} // else
	} // while
	
	return sSearchPath;
} // function