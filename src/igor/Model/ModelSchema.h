// ModelSchema.h ---

#pragma once

/**
 * \file ModelSchema.h
 * \brief The JSON Schemas of the model document, and validation against them.
 *
 * The schemas live in src/igor/Model/schemas/<version>/, one file per definition: the document,
 * its sections, and one file per event type, which the event union refers to. They are embedded
 * in the library at build time and installed under share/igor/schemas/model/. A document names
 * its schema in "$schema", so that an editor validates and completes it.
 *
 * The schema checks each part of a document on its own: keys, types, enumerations, ranges. What
 * relates one part to another (the conditioning naming events, a template for each label, the
 * parameters' axes matching the conditioning, dense indices) is the readers' (ModelJson.h).
 */

#include <igor/Model/Export.h>

#include <nlohmann/json_fwd.hpp>

#include <initializer_list>
#include <string>
#include <string_view>
#include <vector>

namespace igor::model::schema {

/// Where the schemas of this version are published: each schema's "$id" is this plus its file name.
inline constexpr const char *kBaseUri =
        "https://raw.githubusercontent.com/IGoR-immuno-stack/IGoR/master/src/igor/Model/schemas/0.3/";

/// The "$schema" a model document carries: its schema's "$id".
MODEL_EXPORT std::string document_schema_uri();

/// The embedded schema files, by name ("model.schema.json", "event-Deletion.schema.json", ...).
MODEL_EXPORT std::vector<std::string> schema_files();

/// The text of an embedded schema file. \throws std::out_of_range for a name it does not have.
MODEL_EXPORT std::string_view schema_text(const std::string &file);

/**
 * \brief Check \a instance against the schema of \a file, which must also have the keys \a required.
 *
 * \a what names the instance in the message. The validators are built once, on first use.
 * \throws std::runtime_error listing every violation, each at its JSON pointer.
 */
MODEL_EXPORT void validate(const nlohmann::json &instance, const std::string &file, const std::string &what,
                           std::initializer_list<const char *> required = {});

}  // namespace igor::model::schema

//
// ModelSchema.h ends here
